struct GaussianSqrt
    mean
    SqrtCovar
end

rand(g::GaussianSqrt) = g.mean + g.SqrtCovar * randn(size(g.SqrtCovar,2))


struct Gaussian
    mean
    covar
end

GaussianSqrt(g::Gaussian) = GaussianSqrt(g.mean,cholesky(g.covar).U)
rand(g::Gaussian) = rand(GaussianSqrt(g))



struct VectorFun{T,F} <: AbstractVector{T}
    len::Int
    fun::F
end

import Base: size, getindex

Base.size(v::VectorFun) = (v.len,)
Base.getindex(v::VectorFun,n::Integer) = v.fun(n)


VectorFun(T,len,fun) = VectorFun{T,typeof(fun)}(len,fun)
