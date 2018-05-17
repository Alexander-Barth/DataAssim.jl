import DataAssim
using Test
using LinearAlgebra

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
Xf = reshape(sin.(3*(1:(n*N))),n,N)
H = reshape(1:(m*n),m,n)

R = Matrix(2*I,m,m)

xf = mean(Xf, dims=2)
Xfp = Xf .- xf


Pf = (Xfp * Xfp') / (N-1)
K = Pf * H'*inv(H*Pf*H' + R)
Pa_check = Pf - K*H*Pf
xa_check = xf + K*(y - H*xf)

methods = [DataAssim.EnSRF, DataAssim.EAKF, DataAssim.ETKF, DataAssim.ETKF2, DataAssim.SEIK, DataAssim.ESTKF, DataAssim.serialEnSRF]

# Non-serial algorithms (all observations at once) with H

for method in methods
    #Xa,xa = method(Xf,H*Xf,y,R; debug=debug, tolerance=tol, H=H)
    Xa,xa = method(Xf,H*Xf,y,R,H; debug=debug, tolerance=tol)
    Xap = Xa .- view(xa,:,1:1)

    # check analysis
    @test xa ≈ xa_check

    # check analysis ensemble mean
    @test mean(Xa, dims = 2) ≈ xa_check

    # check analysis ensemble variance
    @test (Xap * Xap') / (N-1) ≈ Pa_check
end

# check if no error, output is stochastic
Xa,xa = DataAssim.EnKF(Xf,H*Xf,y,R,H; debug=debug, tolerance=tol)

# local analysis

part = collect(1:n)

selectObs(i) = ones(m)

Xag,xa = DataAssim.ETKF2(Xf,H*Xf,y,R,H)
# local analysis
Xal,xal = DataAssim.local_ETKF2(Xf,H,y,diag(R),part,selectObs)

# check of global is the same as local analysis if all weights are 1
@test Xag ≈ Xal

# test compact function

@test DataAssim.compact_locfun(0) ≈ 1 atol=tol
# reference value from octave code
@test DataAssim.compact_locfun(1.5) ≈ 19/1152 atol=tol
@test DataAssim.compact_locfun(2) ≈ 0 atol=tol
@test DataAssim.compact_locfun(3) ≈ 0 atol=tol
