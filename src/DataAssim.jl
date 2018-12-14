# Functions for data assimilation generated in the frame of the Sangoma project


module DataAssim
using Test
using LinearAlgebra
using Printf
using Statistics
using DIVAnd
import DIVAnd: pack, unpack

export compact_locfun

"""
Abstract base-class of models. A model should implement forecast step, 
tangent-linear and adjoint step
"""
abstract type AbstractModel
end
export AbstractModel

function tgl(M::AbstractModel,t,x,dx::AbstractVecOrMat)
    dx2 = similar(dx)
    for i = 1:size(dx,2)
        dx2[:,i] = tgl(M,t,x,dx[:,i])
    end
    return dx2
end

mutable struct ModelMatrix{T <: Union{AbstractMatrix,UniformScaling}} <: AbstractModel
    M::T
end

(M::ModelMatrix)(t,x,η = zeros(size(x))) = M.M*x + η
tgl(M::ModelMatrix,t,x,dx::AbstractVecOrMat) = M.M*dx
adj(M::ModelMatrix,t,x,dx::AbstractVecOrMat) = M.M'*dx

struct ModelFun{F,F2,F3} <: AbstractModel
    forecast::F
    tgl::F2
    adj::F3
end

(M::ModelFun)(t,x,η = zeros(size(x))) = M.forecast(t,x,η)
tgl(M::ModelFun,t,x,dx::AbstractVector) = M.tgl(t,x,dx)
adj(M::ModelFun,t,x,dx::AbstractVector) = M.adj(t,x,dx)


export tgl, adj

for method = [:EnSRF, :EAKF, :ETKF, :ETKF2, :SEIK, :ESTKF, :serialEnSRF, :EnKF]

    @eval begin
        """
        Xa,xa = "ensemble_analysis"(Xf,HXf,y,R,H,...)
    
    Computes analysis ensemble Xa based on forecast ensemble Xf
    and observations y using various ensemble scheme.
    The function name "ensemble_analysis" can be EnSRF, EAKF, ETKF, SEIK, ESTKF, serialEnSRF or EnKF.

# Input arguments:
* `Xf`: forecast ensemble (n x N)
* `HXf`: the observation operator applied on the ensemble (product H*Xf)
* `y`: observations (m)
* `R`: observation error covariance  (m x m).
* `H`: operator (m x n). Except for the serialEnSRF it is never used and can be empty

# Optional keywords arguments:
* `debug`: set to true to enable debugging. Default (false) is no debugging.
* `tolerance`: expected rounding error (default 1e-10) for debugging checks. This is not used if debug is false.

# Output arguments:
* `Xa`: the analysis ensemble (n x N)
* `xa`: the analysis ensemble mean (n)

Notations follows:
Sangoma D3.1 http://data-assimilation.net/Documents/sangomaDL3.1.pdf
"""       
        function $method(Xf,HXf,y,R,H; debug = false, tolerance=1e-10)
            #function $method(Xf,HXf,y,R,H)
            #debug = false; tolerance=1e-10;
            tol = tolerance
            
            # ensemble size
            N = size(Xf,2)
            
            # number of observations
            m = size(y,1)
            
            xf = mean(Xf,dims = 2)[:,1]
            Xfp = Xf - repeat(xf, inner = (1,N))
            
            Hxf = mean(HXf, dims = 2)[:,1]
            S = HXf - repeat(Hxf, inner = (1,N))

            F = S*S' + (N-1) * R

            if debug
                Pf = (Xfp * Xfp') / (N-1)
                HPfH = (S * S') / (N-1)
                PfH = (Xfp * S') / (N-1)

                K = PfH * inv(HPfH + R)
            end

            if $method == EnSRF
                # EnSRF
                e = eigen(Symmetric(F))
                Gamma_S = e.vectors
                Lambda_S = Diagonal(e.values)
                
                #Lambda_S,Gamma_S = eig(F)
                #Lambda_S = Diagonal(Lambda_S)
                
                X_S = S'*Gamma_S * sqrt.(inv(Lambda_S))
                U_S,Sigma_S,Z_S = svd(X_S)

                Sigma_S = Diagonal(Sigma_S)

                Xap = Xfp * (U_S * (sqrt.(I - Sigma_S*Sigma_S') * U_S'))
                xa = xf + Xfp * (S' * (Gamma_S * (Lambda_S \ (Gamma_S' * (y - Hxf)))))

            elseif $method == serialEnSRF
                # EnSRF with serial observation processing
                for iobs = 1:m
                    # the number 1 with the same element type of Xf
                    one = eltype(Xf)(1)
                    
                    # H[[iobs],:] is necessary instead of H[iobs,:] to make it a row vector
                    Hloc = H[[iobs],:]
                    yloc = y[[iobs]]

                    Sloc = Hloc*Xfp
                    Hxfloc = Hloc*xf

                    # ()[1] makes a scalar instead of a vector of size 1
                    Floc = (Sloc*Sloc')[1] + (N-1)*R[iobs, iobs] 

                    Kloc = Xfp*Sloc' / Floc

                    xa = xf + Kloc * (yloc - Hxfloc)
                    alpha = one / (one + sqrt( (N-1)*R[iobs,iobs]/Floc) )
                    Xap = Xfp - alpha * Kloc * Sloc

                    Xfp = Xap
                    xf = xa
                end
            elseif $method == ETKF
                # ETKF with decomposition of Stilde
                sqrtR = sqrt(R)
                Stilde = sqrt(1/(N-1)) * (sqrtR \ S)

                # "economy size" SVD decomposition
                #U_T,Sigma_T,V_T = svd(Stilde')
                U_T,Sigma_T,V_T = svd(copy(Stilde'))
                Sigma_T = Diagonal(Sigma_T)

                if size(Sigma_T,2) > N
                    Sigma_T = Sigma_T(:,1:N)
                    V_T = V_T(:,1:N)
                end
                Ndim = size(Sigma_T,1)

                TTt = I - S'*(F\S)

                if debug
                    # ETKF-TTt
                    @test TTt ≈ U_T * ((I + Sigma_T*Sigma_T') \ U_T')

                    K2 = 1/sqrt(N-1) * Xfp *
                        (Stilde' * ((Stilde*Stilde'+ I) \ inv(sqrtR)))

                    @test K ≈ K2
                end

                Xap = Xfp * (U_T * (sqrt.(I + Sigma_T*Sigma_T') \ U_T'))
                xa = xf + 1/sqrt(N-1) *  Xfp * (U_T * ((I + Sigma_T'*Sigma_T) \
                                                       (Sigma_T * V_T' * (sqrtR \ (y - Hxf)))))

            elseif $method == ETKF2
                # ETKF with square-root of invTTt (e.g. Hunt et al., 2007)

                invR_S = R \ S
                invTTt = Symmetric((N-1) * I + S' * invR_S)

                # eig is symmetric, thus the type of its eigenvalues are known

                e = eigen(invTTt)
                U_T = e.vectors
                Sigma_T = Diagonal(e.values)


                if debug
                    # ETKF2-eig

                    @test U_T*Sigma_T * U_T' ≈ invTTt
                end

                T = U_T * (sqrt.(Sigma_T) \ U_T')

                if debug
                    # ETKF2-sym. square root
                    @test T*T ≈ inv(invTTt)
                end

                Xap = sqrt(N-1) * Xfp * T
                xa = xf + Xfp * (U_T * (inv(Sigma_T) * U_T' * (invR_S' * (y - Hxf))))


                # elseif $method == "ETKF3"
                #   # ETKF in style of ESTKF
                #   A_T = zeros(N,N)

                #   for j = 1:N
                #     for i = 1:N
                #       if i == j
                #         A_T(i,j) = 1 - 1/N
                #       else
                #         A_T(i,j) = - 1/N
                #       end
                #     end
                #   end

                #   L = Xf*A_T
                #   HL = HXf*(A_T*L)

                #   if debug
                #     Pf2 = 1/(N-1) * L * ((A_T'*A_T) \ L')
                #     @test Pf,Pf2',"ETKF3-Pf",tol)
                #     Pf2 = 1/(N-1) * L * L'
                #     @test Pf,Pf2',"ETKF3-Pf2",tol)
                #   end

                #   invTTt = (N-1)*eye(N) + HL' * (R \ HL)

                #   Sigma_T,U_T = eig(invTTt)
                #   Sigma_T = Diagonal(Sigma_T)

                #   T = U_T * (sqrt(Sigma_T) \ U_T')
                #   #T = sqrt(inv(invTTt))

                #   if debug
                #     @test T*T,inv(invTTt),"ETKF3-sym. square root",tol)
                #   end

                #   Xap = sqrt(N-1) * L*T * A_T'
                #   xa = xf + L * (U_T * (inv(Sigma_T) * U_T' * (HL' * (R \ (y - Hxf)))))

            elseif $method == EAKF
                # EAKF
                sqrtR = sqrt(R)

                Stilde = sqrt(1/(N-1)) * (sqrtR \ S)

                # Sigma_A should have the size (N-1) x (N-1)

# issue
# https://github.com/JuliaLang/julia/issues/27132

U_A,Sigma_A,V_A = svd(Stilde')
Sigma_A = Diagonal(Sigma_A)

if debug
    # last singular should be zero
    @test abs.(Sigma_A[N,N]) ≈ 0 atol=tol
end

U_A = U_A[:,1:N-1]
Sigma_A = Sigma_A[1:N-1,1:N-1]
V_A = V_A[:,1:N-1]

# eigenvalue decomposition of Pf
# Gamma_A should have the size (N-1) x (N-1)

#[Z_A,Gamma_A] = eig(Pf)
Z_A,sqrtGamma_A,dummy = svd(Xfp)
sqrtGamma_A = Diagonal(sqrtGamma_A)

if debug
    # last eigenvalue should be zero
    @test abs(sqrtGamma_A[end,end]) ≈ 0 atol=tol
end

Gamma_A = sqrtGamma_A.^2 / (N-1)
Gamma_A = Gamma_A[1:N-1,1:N-1]
Z_A = Z_A[:,1:N-1]

if debug
    # EAKF-decomposition
    @test Z_A * Gamma_A * Z_A' ≈ Pf
end

Xap = 1/sqrt(N-1) * Xfp * (U_A * ((sqrt.(I + Sigma_A'*Sigma_A)) \
                                  (sqrt.(inv(Gamma_A)) * (Z_A' * Xfp))))

xa = xf + 1/sqrt(N-1) *  Xfp * (U_A * ((I + Sigma_A'*Sigma_A) \
                                       (Sigma_A * V_A' * (sqrtR \ (y - Hxf)))))

elseif $method == SEIK
# SEIK

A = [I; zeros(N-1)']
A = A - ones(N,N-1)/N

L = Xf*A
HL = HXf*A

if debug
    # SEIK-L matrix
    @test L ≈ Xfp[:,1:N-1]

    # SEIK-Pf
    Pf2 = 1/(N-1) * L * ((A'*A) \ L')
    @test Pf ≈ Pf2
end

invTTt = (N-1)*(A'*A) + HL' * (R \ HL)
TTt = inv(invTTt)

T = cholesky(inv(invTTt)).U'

if debug
    # SEIK-Cholesky decomposition
    @test T*T' ≈ inv(invTTt)
end

# add omega
w = ones(N,1)/sqrt(N)
Omega = householder(w)
Xap = sqrt(N-1) * L*T*Omega'

xa = xf + L * (TTt * (HL' * (R \ (y - Hxf))))
elseif $method == ESTKF
# ESTKF
A = zeros(N,N-1)

for j = 1:N-1
    for i = 1:N
        if i == j && i < N
            A[i,j] = 1 - 1/N * 1/(1/sqrt(N)+1)
        elseif i < N
            A[i,j] = - 1/N * 1/(1/sqrt(N)+1)
        else
            A[i,j] = - 1/sqrt(N)
        end
    end
end

L = Xf*A
HL = HXf*A

if debug
    #  ESTKF-Pf
    Pf2 = 1/(N-1) * L * ((A'*A) \ L')
    @test Pf ≈ Pf2

    # ESTKF-Pf2
    Pf2 = 1/(N-1) * L * L'
    @test Pf ≈ Pf2
end

invTTt = (N-1)*I + HL' * (R \ HL)

Sigma_E,U_E = eigen(invTTt)
Sigma_E = Diagonal(Sigma_E)

T = U_E * (sqrt.(Sigma_E) \ U_E')
#T = sqrt(inv(invTTt))

if debug
    # ESTKF-sym. square root
    @test T*T ≈ inv(invTTt)
end

Xap = sqrt(N-1) * L*T * A'
xa = xf + L * (U_E * (inv(Sigma_E) * U_E' * (HL' * (R \ (y - Hxf)))))
elseif $method == EnKF
# EnKF
sqrtR = sqrt(R)

# perturbation of observations
Yp = sqrtR * randn(m,N)
Y = Yp + repeat(y,inner=(1,N))

U_F,Sigma_F,V_F = svd(S + Yp)
Sigma_F = Diagonal(Sigma_F)

G_F,Gamma_F,Z_F = svd(S'*(U_F*inv(Sigma_F)))

# unclear in manuscript
#Xap = Xfp * G_F * sqrt(I - Gamma_F*Gamma_F') * G_F'

Xa = Xf + Xfp * (S' * (U_F * (inv(Sigma_F)^2 * (U_F' * (Y-HXf)))))

xa = mean(Xa, dims = 2)
Xap = Xa - repeat(xa,inner=(1,N))
end

Xa = Xap + repeat(xa,inner=(1,N))

return Xa,xa

end # function $method

export $method

"""
[Xa,xa] = local_ensemble_analysis(...
   Xf,H,y,diagR,part,selectObs,method,...)

Computes analysis ensemble Xa based on forecast ensemble Xf using the
observation y.

Inputs:
Xf: forecast ensemble (n x N)
H: observation operator (m x n)
y: observation (m x 1)
diagR: diagonal of the observation error covariance R (m x 1)
part: vector of integer "labels". Every element of the state vector with the
  same number belong to the same subdomain
selectObs: callback routine to select observations with a within a subdomain.
  As input is takes an integer representing the index of the state vector and
  returns a vector of weights (m x 1).
  For example:
     selectObs(i) = exp(- ((x[i] - xobs[:]).^2 + (y(i) - yobs[:]).^2)/L^2 );
  or
     selectObs(i) = sangoma_compact_locfun(L,...
         sqrt((x[i] - xobs[:]).^2 + (y[i] - yobs[:]).^2));

  where:
     x and y is the horizontal model grid
     xobs and yobs are the localtion of the observations
     L is a correlation length-scale

method: method is one analysis schemes implemented sangoma_ensemble_analysis
  (except for EnSRF)

Optional inputs:
* display: if true, then display progress (false is the default)
* minweight: analysis is performed using observations for which
   weights is larger than minweight. (default 1e-8)
* HXf: if non empty, then it is the product H Xf. In this case, H is not
   used

Output:
`Xa`: the analysis ensemble (n x N)
`xa`: the analysis ensemble mean (n x 1)

See also: compact_locfun
"""
function $(Symbol("local_" * string(method)))(Xf,H,y,diagR,part,selectObs;
                                              display = false,
                                              minweight = 1e-8,
                                              HXf = [])

    # unique element of partition vector
    p = unique(part);

    Xa = zeros(size(Xf));
    xa = zeros(size(Xf,1),1);

    # do not use isempty here because m might be zero
    if isequal(HXf,[])
        HXf = H*Xf;
    end

    # loop over all zones
    for i=1:length(p)
        if display
            @printf("zone %d out of %d\n",i,length(p));
        end

        sel = findall(part .== p[i]);
        weight = selectObs(sel[1]);

        # restrict to local observations where weight exceeds minweight
        loc = findall(weight .> minweight);
        HXfloc = HXf[loc,:];
        Rloc = Diagonal(diagR[loc] ./ weight[loc]);
        yloc = y[loc];

        Xa[sel,:],xa[sel] = $method(Xf[sel,:],HXfloc,
                                    yloc,Rloc, []);


    end

    return Xa,xa
end # function local_$method

export $(Symbol("local_" * string(method)))

end # @eval begin
end # for method ...


function householder(w)
    n = length(w)-1

    H = zeros(n+1,n)

    w2 = copy(w)
    w2[n+1] = w2[n+1] + sign(w[n+1])

    H = [I; zeros(1,n)] - 1/(abs(w[n+1])+1) *  w2 * w[1:n]'

    return H
end

"""
fun = compact_locfun(r)

Smooth compact localization function described in Gaspari et al. (1999),
equation 4.10.

Input:
  r: matrix of distance (scaled)
Output:
  fun: array of weights (same size as r). fun is zero if r > 2

Reference:
@article {QJ:QJ49712555417,
author = {Gaspari, Gregory and Cohn, Stephen E.},
title = {Construction of correlation functions in two and three dimensions},
journal = {Quarterly Journal of the Royal Meteorological Society},
volume = {125},
number = {554},
publisher = {John Wiley & Sons, Ltd},
issn = {1477-870X},
url = {http://dx.doi.org/10.1002/qj.49712555417},
doi = {10.1002/qj.49712555417},
pages = {723--757},
keywords = {Compactly supported, Convolution, Correlation functions, Data assimilation, Space-limited},
year = {1999},
}
"""
function compact_locfun(r)

    if r <= 1
        return (((-r/4. + 1/2) * r + 5/8) * r - 5/3) * r^2 + 1
    elseif 1 < r <= 2;
        return ((((r/12. - 1/2) * r + 5/8) * r + 5/3) * r - 5) * r + 4 - 2/(3*r);
    else
        return zero(r)
    end
end




include("fourDVar.jl")
include("KalmanFilter.jl")
include("TwinExperiment.jl")
include("lorenz63model.jl")
include("shallow_water1D_model.jl")

export FreeRun, fourDVar, TwinExperiment, LinShallowWater1DModel, KalmanFilter
export AbstractModel, ModelFun
export pack, unpack
export Lorenz63Model
export ModelMatrix
end
