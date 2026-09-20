import DataAssim
using Test
using Statistics
using LinearAlgebra


permutations(t::NTuple{0}) = [t]

# all permutation of tuple t
function permutations(t::T) where T <: Tuple
    result = T[]

    for i in eachindex(t)
        x = t[i]
        rest = (t[1:i-1]..., t[i+1:end]...)

        for p in permutations(rest)
            push!(result, (x, p...))
        end
    end

    return result
end


# n: number of elements in the state vector
# N: ensemble size
# m: number of observations

for (n,N,m) in [permutations((3,6,9))...,(6,6,6)]
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

    method = DataAssim.ETKF

    Xam,xam = method(Xf,H*Xf,y,R,H; debug=debug, tolerance=tol)
    Xap = Xam .- view(xam,:,1:1)

    # check analysis
    @test xam ≈ xa_check

    # check analysis ensemble mean
    @test mean(Xam, dims = 2) ≈ xa_check

    # check analysis ensemble variance
    @test (Xap * Xap') / (N-1) ≈ Pa_check
end
