
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

"""
    ℳ = ModelMatrix(M)

Linear model defined by the matrix `M`.
"""
mutable struct ModelMatrix{T <: Union{AbstractMatrix,UniformScaling}} <: AbstractModel
    M::T
end

(M::ModelMatrix)(t,x,η = zeros(size(x))) = M.M*x + η
tgl(M::ModelMatrix,t,x,dx::AbstractVecOrMat) = M.M*dx
adj(M::ModelMatrix,t,x,dx::AbstractVecOrMat) = M.M'*dx

"""
    ℳ = ModelFun(nonlinear_forecast,tangent_linear_model,adjoint_model)

Model defined by the functions `nonlinear_forecast`,`tangent_linear_model` and 
`adjoint_model`.
"""
struct ModelFun{F,F2,F3} <: AbstractModel
    forecast::F
    tgl::F2
    adj::F3
end

(M::ModelFun)(t,x,η = zeros(size(x))) = M.forecast(t,x,η)
tgl(M::ModelFun,t,x,dx::AbstractVector) = M.tgl(t,x,dx)
adj(M::ModelFun,t,x,dx::AbstractVector) = M.adj(t,x,dx)


export tgl, adj

include("lorenz63model.jl")
include("shallow_water1D_model.jl")
