using DataAssim
using Base.Test

# number of elements in the state vector
n = 10
# ensemble size
N = 3
# number of observations
m = 5

# if debug is true, then internal checks are activated
debug = true

# tolerance for internal checking
tol = 1e-10

# some random data
y = randn(m,1)
Xf = randn(n,N)

H = randn(m,n)

y = 1:m
Xf = reshape(sin(3*(1:(n*N))),n,N)
H = reshape(1:(m*n),m,n)

R = 2*eye(m,m)

xf = mean(Xf,2)
Xfp = Xf - repmat(xf,1,N)


Pf = (Xfp * Xfp') / (N-1)
K = Pf * H'*inv(H*Pf*H' + R)
Pa_check = Pf - K*H*Pf
xa_check = xf + K*(y - H*xf)

method = ["EnSRF","EAKF","ETKF","ETKF2","SEIK","ESTKF","serialEnSRF"]

# Non-serial algorithms (all observations at once) with H

for i = 1:length(method)
  Xa,xa = ensemble_analysis(Xf,H,y,R,method[i],
                                    debug=debug, tolerance=tol)
  Xap = Xa - repmat(xa,1,N)

  # check analysis

  @test xa ≈ xa_check

  # check analysis ensemble mean

  @test mean(Xa,2) ≈ xa_check

  # check analysis ensemble variance

  @test (Xap * Xap') / (N-1) ≈ Pa_check
end

method = ["EnSRF","EAKF","ETKF","ETKF2","SEIK","ESTKF"]

# Non-serial algorithms (all observations at once) with HXf

for i = 1:length(method)
  Xa,xa = ensemble_analysis(Xf,[],y,R,method[i],
                                    debug=debug, tolerance=tol,HXf=H*Xf)
  Xap = Xa - repmat(xa,1,N)

  # check analysis

  @test xa ≈ xa_check

  # check analysis ensemble mean

  @test mean(Xa,2) ≈ xa_check

  # check analysis ensemble variance

  @test (Xap * Xap') / (N-1) ≈ Pa_check
end


# local analysis

# method to use for the analysis
method = "ETKF2";
part = collect(1:n);

selectObs(i) = ones(m);

Xag,xa = ensemble_analysis(Xf,H,y,R,method);
# local analysis
Xal,xal = local_ensemble_analysis(Xf,H,y,diag(R),part,selectObs,method);

# check of global is the same as local analysis if all weights are 1
@test Xag ≈ Xal

# test compact function

@test_approx_eq_eps compact_locfun(0) 1 tol
@test_approx_eq_eps compact_locfun(2) 0 tol
@test_approx_eq_eps compact_locfun(3) 0 tol

