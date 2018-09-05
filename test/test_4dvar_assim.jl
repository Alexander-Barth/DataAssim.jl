using LinearAlgebra
using Test
using Random
using DIVAnd

abstract type AbstractModel
end

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


include("fourDVar.jl")
include("KalmanFilter.jl")
include("TwinExperiment.jl")

include("lorenz63model.jl")
include("shallow_water1D_model.jl")

Random.seed!(12343)


include("test_shallow_water1D_model.jl")

ℳ = Lorenz63Model(0.01)

@test ℳ(0,[1.,2.,3.]) ≈ [1.1065,  2.241665,  2.9430075] atol=1e-3

x = randn(3,10000)
for k = 1:size(x,2)-1
    x[:,k+1] = ℳ(k,x[:,k])
end

check(ℳ,3)

x = randn(4)
ℳ = ModelMatrix(2*I)
@test ℳ(0,x) ≈ 2*x
@test tgl(ℳ,0,0,x) ≈ 2*x
@test adj(ℳ,0,0,x) ≈ 2*x



ℳ = ModelFun((t,x,η) -> 2*x,(t,x,dx) -> 2*dx,(t,x,dx) -> 2*dx)
@test ℳ(0,x) ≈ 2*x
@test tgl(ℳ,0,0,x) ≈ 2*x
@test adj(ℳ,0,0,x) ≈ 2*x



check(ℳ,4)



# test: one obs at IC

n = 2;
m = 1;

#M = @(x) x;
#MT = @(x) x;

M = I
model_fun(t,x,η) = M*x;
model_tgl(t,x,dx) = M*dx;
model_adj(t,x,dx) = M'*dx;
ℳ = ModelMatrix(I)


H = [1 0];
𝓗 = ModelMatrix(H)
xi = [1; 1];
Pi = Matrix(I,n,n)
R = Matrix(I,m,m)

nmax = 0;
yo = randn(m,nmax+1);

# at which time step to assimilate
# 1 is IC, 2 -> after first time step
no=[1];

xa, = fourDVar(xi,Pi,ℳ,yo,R,H,nmax,no);
@inferred fourDVar(xi,Pi,ℳ,yo,R,H,nmax,no)

#[xa3] = pcg(fun,b,xi);

P = Pi;
K = P*H'*inv(H*P*H' + R);
Pa = P - K*H*P;
xa2  = xi + K * (yo - H*xi);

# should be ~0
#rms(test_4dvar_grad(xi,xa2,Pi,M,yo,R,H,nmax,no),zeros(n,1))

# should be ~0
@test xa ≈ xa2

xa3, = KalmanFilter(xi,Pi,ℳ,zeros(size(Pi)),yo,R,H,nmax,no);
# should be ~0
@test xa ≈ xa3

@inferred KalmanFilter(xi,Pi,ℳ,zeros(size(Pi)),yo,R,H,nmax,no)

#-----------------------------------------
# test: two obs at IC (no evolution)

nmax = 1;
yo = randn(m,nmax+1);
yo = [3 7];
no = [1,2];

xa, = fourDVar(xi,Pi,ℳ,yo,R,H,nmax,no);

P = Pi;
K = P*H'*inv(H*P*H' + R);
P = P - K*H*P;
xa2  = xi + K * (yo[:,1] - H*xi);

K = P*H'*inv(H*P*H' + R);
xa2  = xa2 + K * (yo[:,2] - H*xa2);

# should be ~0
@test xa ≈ xa2 atol=1e-14

#𝓗
#𝓜
xa3, = KalmanFilter(xi,Pi,ℳ,zeros(size(Pi)),yo,R,H,nmax,no);
# should be ~0
@test M*xa ≈ xa3[:,end]  atol=1e-14


#-----------------------------------------
# test: one obs at IC, one at next time step (with evolution)

M = [1 -.1; 0.1 1];
model_fun(t,x,η) = M*x;
model_tgl(t,x,dx) = M*dx;
model_adj(t,x,dx) = M'*dx;
ℳ = ModelMatrix(M)

xa, = fourDVar(xi,Pi,ℳ,yo,R,H,nmax,no);
xa2, = KalmanFilter(xi,Pi,ℳ,zeros(size(Pi)),yo,R,H,nmax,no);
# should be ~0
@test M*xa ≈ xa2[:,end] atol=1e-14



#-----------------------------------------
# test: one obs next time step 2 and one at 5
no = [2,5];
nmax = 10;

xa, = fourDVar(xi,Pi,ℳ,yo,R,H,nmax,no);
xa2, = KalmanFilter(xi,Pi,ℳ,zeros(size(Pi)),yo,R,H,nmax,no);
# should be ~0
@test M^(nmax)*xa ≈ xa2[:,end] atol=1e-14


# twin experiment

# test: one obs at IC

n = 2;
m = 1;

M = I
model_fun(t,x,η) = M*x;
model_tgl(t,x,dx) = M*dx;
model_adj(t,x,dx) = M'*dx;
ℳ = ModelMatrix(I)

H = [1 0];
xit = [1; 1];
Pi = Matrix(I,n,n)
R = Matrix(I,m,m)
Q = zeros(n,n);

nmax = 100;

# at which time step to assimilate
# 1 is IC, 2 -> after first time step
no=3:nmax;
method = "4DVar";

xt,xfree,xa,yt,yo = TwinExperiment(ℳ,xit,Pi,Q,R,H,nmax,no,method);

@inferred FreeRun(ℳ,xi,Q,H,nmax,no)
@inferred TwinExperiment(ℳ,xit,Pi,Q,R,H,nmax,no,method)


# lorenz63

ℳ = Lorenz63Model(0.01)

nmax = 20;
no = 5:nmax;
n = 3;
sigma=10;
beta = 8/3;
rho = 28;
dt = 0.02;

xit = [5.; 0.; 0.];
H = [1 0 0];
Pi = Matrix(3*I,n,n)
Q = zeros(n,n);

# $$$ model_fun = @(t,x) lorenz63(x,dt);
# $$$ model_tgl = @(t,x,dx) rungekutta4(0,dx,dt,@(t,dx) [  -sigma, sigma,      0; rho-x(3),    -1,  -x(1); x(2),  x(1),  -beta] * dx);
# $$$ model_adj = @(t,x,dx) rungekutta4(0,dx,dt,@(t,dx) [  -sigma, sigma,      0; rho-x(3),    -1,  -x(1); x(2),  x(1),  -beta]' * dx);
# $$$
# $$$ check_tgl_adj(model,3,0);

method = "4DVar";

#xt,xfree,xa,yt,yo,diag = TwinExperiment(ℳ,xit,Pi,Q,R,H,nmax,no,method);

ℳ = Lorenz63Model(0.05)


if true
nmax = 10000;
#xt,xfree,xa,yt,yo,diag = TwinExperiment(ℳ,xit,Pi,Q,R,H,nmax,no,method);

# true run
xt,yt = FreeRun(ℳ,xit,Q,H,nmax,no);


end

nmax = 100;
no = 5:nmax;
method = "KF";
xt,xfree,xa,yt,yo,diag = TwinExperiment(ℳ,xit,Pi,Q,R,H,nmax,no,method);

if false
    using PyPlot
    subplot(2,1,1)
    plot(xt[1,:],"b",label = "true")
    plot(xfree[1,:],"r",label = "free")
    plot(xa[1,:],"g", label = "assim")
    legend()
    subplot(2,1,2)
    plot(diag.J)
end

