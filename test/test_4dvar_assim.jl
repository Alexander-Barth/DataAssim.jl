using LinearAlgebra
using Test
using Random
using DIVAnd

include("fourDVar.jl")
include("KalmanFilter.jl")
include("TwinExperiment.jl")

Random.seed!(12343)

# test: one obs at IC

n = 2;
m = 1;

#M = @(x) x;
#MT = @(x) x;

M = I
model_fun(t,x,η) = M*x;
model_tgl(t,x,dx) = M*dx;
model_adj(t,x,dx) = M'*dx;

H = [1 0];
xi = [1; 1];
Pi = Matrix(I,n,n)
R = Matrix(I,m,m)

nmax = 0;
yo = randn(m,nmax+1);

# at which time step to assimilate
# 1 is IC, 2 -> after first time step
no=[1];

xa, = fourDVar(xi,Pi,model_fun,model_tgl,model_adj,yo,R,H,nmax,no);
@inferred fourDVar(xi,Pi,model_fun,model_tgl,model_adj,yo,R,H,nmax,no)

#[xa3] = pcg(fun,b,xi);

P = Pi;
K = P*H'*inv(H*P*H' + R);
Pa = P - K*H*P;
xa2  = xi + K * (yo - H*xi);

# should be ~0
#rms(test_4dvar_grad(xi,xa2,Pi,M,yo,R,H,nmax,no),zeros(n,1))

# should be ~0
@test xa ≈ xa2

xa3, = KalmanFilter(xi,Pi,model_fun,model_tgl,zeros(size(Pi)),yo,R,H,nmax,no);
# should be ~0
@test xa ≈ xa3

@inferred KalmanFilter(xi,Pi,model_fun,model_tgl,zeros(size(Pi)),yo,R,H,nmax,no)

#-----------------------------------------
# test: two obs at IC (no evolution)

nmax = 1;
yo = randn(m,nmax+1);
yo = [3 7];
no = [1,2];

xa, = fourDVar(xi,Pi,model_fun,model_tgl,model_adj,yo,R,H,nmax,no);

P = Pi;
K = P*H'*inv(H*P*H' + R);
P = P - K*H*P;
xa2  = xi + K * (yo[:,1] - H*xi);

K = P*H'*inv(H*P*H' + R);
xa2  = xa2 + K * (yo[:,2] - H*xa2);

# should be ~0
@test xa ≈ xa2 atol=1e-14


xa3, = KalmanFilter(xi,Pi,model_fun,model_tgl,zeros(size(Pi)),yo,R,H,nmax,no);
# should be ~0
@test M*xa ≈ xa3[:,end]  atol=1e-14


#-----------------------------------------
# test: one obs at IC, one at next time step (with evolution)

M = [1 -.1; 0.1 1];
model_fun(t,x,η) = M*x;
model_tgl(t,x,dx) = M*dx;
model_adj(t,x,dx) = M'*dx;

xa, = fourDVar(xi,Pi,model_fun,model_tgl,model_adj,yo,R,H,nmax,no);
xa2, = KalmanFilter(xi,Pi,model_fun,model_tgl,zeros(size(Pi)),yo,R,H,nmax,no);
# should be ~0
@test M*xa ≈ xa2[:,end] atol=1e-14



#-----------------------------------------
# test: one obs next time step 2 and one at 5
no = [2,5];
nmax = 10;

xa, = fourDVar(xi,Pi,model_fun,model_tgl,model_adj,yo,R,H,nmax,no);
xa2, = KalmanFilter(xi,Pi,model_fun,model_tgl,zeros(size(Pi)),yo,R,H,nmax,no);
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

xt,xfree,xa,yt,yo = TwinExperiment(model_fun,model_tgl,model_adj,xit,Pi,Q,R,H,nmax,no,method);

@inferred FreeRun(model_fun,xi,Q,H,nmax,no)

#=

# lorenz63

nmax = 20;
no = 5:nmax;
n = 3;
sigma=10;
beta = 8/3;
rho = 28;
dt = 0.02;

xit = [5 0 0]';
H = [1 0 0];
Pi = 3*eye(n);
Q = zeros(n,n);

# $$$ model_fun = @(t,x) lorenz63(x,dt);
# $$$ model_tgl = @(t,x,dx) rungekutta4(0,dx,dt,@(t,dx) [  -sigma, sigma,      0; rho-x(3),    -1,  -x(1); x(2),  x(1),  -beta] * dx);
# $$$ model_adj = @(t,x,dx) rungekutta4(0,dx,dt,@(t,dx) [  -sigma, sigma,      0; rho-x(3),    -1,  -x(1); x(2),  x(1),  -beta]' * dx);
# $$$
# $$$ check_tgl_adj(model,3,0);

method = '4DVar';

#[xt,xfree,xa,yt,yo,diag] = TwinExperiment(model_fun,model_tgl,model_adj,xit,Pi,Q,R,H,nmax,no,method);


if 0
model_fun = @(t,x) x + [sigma*(x(2)-x(1)); x(1)*(rho-x(3)) - x(2); x(1)*x(2) - beta * x(3)]*dt;
model_tgl = @(t,x,dx) dx + [  -sigma, sigma,      0; rho-x(3),    -1,  -x(1); x(2),  x(1),  -beta] * dx*dt;
model_adj = @(t,x,dx) dx + [  -sigma, sigma,      0; rho-x(3),    -1,  -x(1); x(2),  x(1),  -beta]' * dx*dt;
check_tgl_adj(model,3,0);
end

model = lorenz63model(0.05);
check_tgl_adj(model,3,0);

if 1
nmax = 10000;
#[xt,xfree,xa,yt,yo,diag] = TwinExperiment(model_fun,model_tgl,model_adj,xit,Pi,Q,R,H,nmax,no,method);

# true run
[xt,yt] = FreeRun(model_fun,xit,Q,H,nmax,no);
rg(xt)

end
nmax = 100;
no = 5:nmax;
method = 'KF';
[xt,xfree,xa,yt,yo,diag] = TwinExperiment(model_fun,model_tgl,model_adj,xit,Pi,Q,R,H,nmax,no,method);

if 0
subplot(2,1,1)
hold on
plot(xt[1,:]','b')
plot(xfree[1,:]','r')
plot(xa[1,:]','g')
legend('true','free','assim');
subplot(2,1,2)
plot(diag.J)
end
=#
