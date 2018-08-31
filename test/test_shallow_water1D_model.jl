

dt = 1.
g = 9.81
h = 100.
imax = 101
L = 10000
ℳ = LinShallowWater1DModel(dt,g,h,L,imax)


n = 2*imax-1
check(ℳ,n)
