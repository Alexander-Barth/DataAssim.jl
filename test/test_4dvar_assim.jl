using LinearAlgebra
using Test
using Random
#using DIVAnd
using DataAssim

function check(ℳ::AbstractModel,n,t = 0,ϵ = 1e-5)
    dx = randn(n)
    x = randn(n)
    dx2 = randn(n)

    @test (ℳ(t,x + ϵ*dx) - ℳ(t,x - ϵ*dx)) / (2*ϵ)  ≈ tgl(ℳ,t,x,dx) atol=10*ϵ^2
    @test dx2 ⋅ tgl(ℳ,t,x,dx) ≈ adj(ℳ,t,x,dx2) ⋅ dx   atol=1e-7

    dX = randn(n,3)
    MdX = tgl(ℳ,t,x,dX)
    @test tgl(ℳ,t,x,dX[:,1]) ≈ MdX[:,1]
end

Random.seed!(12343)

#include("test_shallow_water1D_model.jl")

@testset "Lorenz63Model ajoint" begin
    ℳ = Lorenz63Model(0.01)
    @test ℳ(0,[1.,2.,3.]) ≈ [1.1065,  2.241665,  2.9430075] atol=1e-3
    check(ℳ,3)
end

@testset "model matrix" begin
    x = randn(4)
    ℳ = ModelMatrix(2*I)
    @test ℳ(0,x) ≈ 2*x
    @test tgl(ℳ,0,0,x) ≈ 2*x
    @test adj(ℳ,0,0,x) ≈ 2*x
end

@testset "model function" begin
    x = randn(4)
    ℳ = ModelFun((t,x,η) -> 2*x,(t,x,dx) -> 2*dx,(t,x,dx) -> 2*dx)
    @test ℳ(0,x) ≈ 2*x
    @test tgl(ℳ,0,0,x) ≈ 2*x
    @test adj(ℳ,0,0,x) ≈ 2*x
    check(ℳ,4)
end

@testset "4DVar (one observation at IC)" begin
    n = 2
    m = 1
    M = I
    ℳ = ModelMatrix(M)
    H = [1 0]
    𝓗 = ModelMatrix(H)
    xi = [1; 1]
    Pi = Matrix(I,n,n)

    nmax = 0
    yo = [randn(m) for i = 1:nmax+1]
    R = fill(Matrix(I,m,m),nmax+1)

    # at which time step to assimilate
    # 1 is IC, 2 -> after first time step
    no = [1]

    xa, = fourDVar(xi,Pi,ℳ,yo,R,𝓗,nmax,no)

    P = Pi
    K = P*H'*inv(H*P*H' + R[1])
    Pa = P - K*H*P
    xa2  = xi + K * (yo[1] - H*xi)
    @test xa ≈ xa2

    xa3, = KalmanFilter(xi,Pi,ℳ,zeros(size(Pi)),yo,R,𝓗,nmax,no)
    @test xa ≈ xa3
end


@testset "4DVar (two obs at IC; no evolution)" begin
    n = 2
    m = 1
    xi = [1; 1]
    Pi = Matrix(I,n,n)
    M = I
    ℳ = ModelMatrix(M)
    H = [1 0]
    𝓗 = ModelMatrix(H)


    nmax = 1
    yo = [[3], [7]]
    R = fill(Matrix(I,m,m),nmax+1)
    no = [1,2]

    xa, = fourDVar(xi,Pi,ℳ,yo,R,𝓗,nmax,no)

    P = Pi
    K = P*H'*inv(H*P*H' + R[1])
    P = P - K*H*P
    xa2  = xi + K * (yo[1] - H*xi)

    K = P*H'*inv(H*P*H' + R[2])
    xa2  = xa2 + K * (yo[2] - H*xa2)

    @test xa ≈ xa2 atol=1e-14

    xa3, = KalmanFilter(xi,Pi,ℳ,zeros(size(Pi)),yo,R,𝓗,nmax,no)
    @test M*xa ≈ xa3[:,end]  atol=1e-14
end

#-----------------------------------------
# test: one obs at IC, one at next time step (with evolution)

@testset "4DVar with evolution" begin
    n = 2
    m = 1
    H = [1 0]
    𝓗 = ModelMatrix(H)
    xi = [1; 1]
    Pi = Matrix(I,n,n)
    nmax = 1
    no = [1]

    yo = [randn(m) for i = 1:nmax+1]
    R = fill(Matrix(I,m,m),nmax+1)

    M = [1 -.1; 0.1 1]
    ℳ = ModelMatrix(M)

    xa, = fourDVar(xi,Pi,ℳ,yo,R,𝓗,nmax,no)
    xa2, = KalmanFilter(xi,Pi,ℳ,zeros(size(Pi)),yo,R,𝓗,nmax,no)
    # should be ~0
    @test M*xa ≈ xa2[:,end] atol=1e-10

end





#-----------------------------------------
# test: one obs next time step 2 and one at 5
@testset "4DVar4" begin
    n = 2
    m = 1
    H = [1 0]
    𝓗 = ModelMatrix(H)
    xi = [1; 1]
    Pi = Matrix(I,n,n)
    R = i -> Matrix(I,m,m)
    no = [2,5]

    yo = [randn(m) for i = 1:length(no)]
    R = fill(Matrix(I,m,m),length(no))

    nmax = 10
    M = [1 -.1; 0.1 1]
    ℳ = ModelMatrix(M)

    xa, = fourDVar(xi,Pi,ℳ,yo,R,𝓗,nmax,no)
    xa2, = KalmanFilter(xi,Pi,ℳ,zeros(size(Pi)),yo,R,𝓗,nmax,no)
    # should be ~0
    @test M^(nmax)*xa ≈ xa2[:,end] atol=1e-10
end

#-----------------------------------------
# test: one obs next time step 2 and one at 5
@testset "KalmanFilter Function API" begin
    n = 2
    m = 1
    H = [1 0]
    𝓗 = ModelMatrix(H)
    xi = [1; 1]
    Pi = Matrix(I,n,n)
    no = [2,5]

    nmax = 10
    M = [1 -.1; 0.1 1]
    ℳ = ModelMatrix(M)

    yo = [randn(m) for i = 1:length(no)]
    R = fill(Matrix(I,m,m),length(no))

    xa, = fourDVar(xi,Pi,ℳ,yo,R,𝓗,nmax,no)
    xa2, = KalmanFilter(xi,Pi,ℳ,zeros(size(Pi)),
                        yo,R,𝓗,nmax,no)
    @test M^(nmax)*xa ≈ xa2[:,end] atol=1e-10
end

@testset "twin experiment (no evolution)" begin
    xi = [1.; 1.]
    n = 2
    m = 1

    M = I
    ℳ = ModelMatrix(I)

    H = [1 0]
    𝓗 = ModelMatrix(H)
    xit = [1.; 1.]
    Pi = Matrix(I,n,n)
    R = i -> Matrix(I,m,m)
    Q = zeros(n,n)

    nmax = 100

    # at which time step to assimilate
    # 1 is IC, 2 -> after first time step
    no = 3:nmax
    R = fill(Matrix(I,m,m),length(no))

    method = "4DVar"

    xt,xfree,xa,yt,yo = TwinExperiment(ℳ,xit,Pi,Q,R,𝓗,nmax,no,method)

    @inferred FreeRun(ℳ,xi,Q,𝓗,nmax,no)
    @inferred TwinExperiment(ℳ,xit,Pi,Q,R,𝓗,nmax,no,method)

    @test_throws ErrorException TwinExperiment(ℳ,xit,Pi,Q,R,𝓗,nmax,no,"best method")

end



@testset "twin experiment: Lorenz 63 " begin
    m = 1

    ℳ = Lorenz63Model(0.01)

    nmax = 20
    no = 5:nmax
    n = 3
    R = fill(Matrix(I,m,m),length(no))


    xit = [5.; 0.; 0.]
    H = [1 0 0]
    𝓗 = ModelMatrix(H)
    Pi = Matrix(3*I,n,n)
    Q = zeros(n,n)

    method = "4DVar"
    xt,xfree,xa,yt,yo,diag_ = TwinExperiment(ℳ,xit,Pi,Q,R,𝓗,nmax,no,method)
    @test norm(xt[:,end] - xa[:,end]) ≈ 0 atol=3

    ℳ = Lorenz63Model(0.05)

    if true
        nmax = 10000
        #xt,xfree,xa,yt,yo,diag_ = TwinExperiment(ℳ,xit,Pi,Q,R,H,nmax,no,method)

        # true run
        xt,yt = FreeRun(ℳ,xit,Q,𝓗,nmax,no)
    end

    nmax = 100
    no = 5:nmax
    R = fill(Matrix(I,m,m),length(no))
    method = "KF"
    xt,xfree,xa,yt,yo,diag_ = TwinExperiment(ℳ,xit,Pi,Q,R,𝓗,nmax,no,method)

    @test norm(xt[:,end] - xa[:,end]) ≈ 0 atol=3

    if false
        using PyPlot
        plot(xt[1,:],"b",label = "true")
        plot(xfree[1,:],"r",label = "free")
        plot(xa[1,:],"g", label = "assim")
        legend()
        title("Extended Kalman filter applied to the Lorenz-63 model")
        savefig("EKF-Lorenz63.svg")
    end
end


