using SparseArrays

dt = 1.
g = 9.81
h = 100.
imax = 101
L = 10000
ℳ = LinShallowWater1DModel(dt,g,h,L,imax)



n = 2*imax-1
check(ℳ,n)

ζt = exp.( - (ℳ.x_r .- L/2).^2 )
ut = zeros(size(ℳ.x_u))

xit = pack(ζt,ut)
Pi = Diagonal(ones(n))

Q = Diagonal(zeros(n))
R = Diagonal(ones(size(ζt)))
method = "4DVar";
H = sparse(1:imax-1,1:imax-1,ones(imax-1),imax-1,n)
nmax = 100

no=3:5:nmax;

xt,xfree,xa,yt,yo = TwinExperiment(ℳ,xit,Pi,Q,R,H,nmax,no,method);

