#=
 +-------------+-------------+-------------+-------------+-------------+
 |      *      |      *      |      *      |      *      |      *      |
 +-------------+-------------+-------------+-------------+-------------+
u[1]   ζ[1]   u[2]   ζ[2]   u[3]    .   .   .   .   .   .   ζ[imax-1] u[imax]

=#


mutable struct LinShallowWater1DModel{T} <: AbstractModel
    dt::T
    g::T
    h::T
    x_r::Vector{T}
    x_u::Vector{T}
end

"""

## Example

```julia
dt = 1.
g = 9.81
h = 100
imax = 101
L = 10000
LinShallowWater1DModel(dt,g,h,L,imax)
```
"""
function LinShallowWater1DModel(dt,g,h,L,imax)
    x_u = collect(range(0,stop=L,length=imax))
    x_r = (x_u[2:end] + x_u[1:end-1])/2
    return LinShallowWater1DModel(dt,g,h,x_r,x_u)
end

function pack(ζ,u)
    x = zeros(eltype(ζ),(length(ζ)+length(u)))
    x[1:length(ζ)] = ζ
    x[length(ζ)+1:end] = u
    return x
end

function unpack(x)
    imax = (length(x)-1) ÷ 2
    return x[1:imax],x[imax+1:end]
end

function (ℳ::LinShallowWater1DModel)(t,x,eta = zeros(eltype(x),size(x)))
    ζ,u = unpack(x)

    for i=1:length(ζ)
      ζ[i] = ζ[i] - ℳ.dt * ℳ.h * (u[i+1] - u[i])/((ℳ.x_u[i+1] - ℳ.x_u[i]))
    end

    for i = 2:length(u)-1
      u[i] = u[i] - ℳ.dt * ℳ.g * (ζ[i] - ζ[i-1])/((ℳ.x_r[i] - ℳ.x_r[i-1]))
    end

    return pack(ζ,u)
end

tgl(ℳ::LinShallowWater1DModel,t,x,dx::AbstractVector) = ℳ(t,dx,dx)

function adj(ℳ::LinShallowWater1DModel,t,x,dx::AbstractVector)
    ζ,u = unpack(dx)

    # first transform the momentum equation
    for i = 2:length(u)-1
        #u[i] = u[i] - ℳ.dt * g * (ζ[i] - ζ[i-1])/((ℳ.x_r[i] - ℳ.x_r[i-1]))
        ζ[i] += - ℳ.dt * ℳ.g * u[i]/((ℳ.x_r[i] - ℳ.x_r[i-1]))
        ζ[i-1] +=  ℳ.dt * ℳ.g * u[i]/((ℳ.x_r[i] - ℳ.x_r[i-1]))
    end

    for i=1:length(ζ)
        u[i+1] += - ℳ.dt * ℳ.h *  ζ[i] /((ℳ.x_u[i+1] - ℳ.x_u[i]))
        u[i] +=  ℳ.dt * ℳ.h *  ζ[i] /((ℳ.x_u[i+1] - ℳ.x_u[i]))
    end

    return pack(ζ,u)
end

