mutable struct Lorenz63Model{T} <: AbstractModel
    dt::T
    σ::T
    β::T
    ρ::T
end


Lorenz63Model(dt,σ=10.,β = 8/3.,ρ = 28.) = Lorenz63Model(dt,σ,β,ρ)

function (ℳ::Lorenz63Model)(t,x,eta = zeros(eltype(x),size(x)))
    return rungekutta2(t,x,ℳ.dt,(t,x) -> lorenz63_dxdt(ℳ.σ,ℳ.β,ℳ.ρ,x)) + eta
end

function tgl(ℳ::Lorenz63Model,t,x,dx::AbstractVector)
    f_tgl(t,x,dx)  = lorenz63_dxdt_tgl(ℳ.σ,ℳ.β,ℳ.ρ,x,dx)
    f(t,x) = lorenz63_dxdt(ℳ.σ,ℳ.β,ℳ.ρ,x)
    return rungekutta2_tgl(t,x,ℳ.dt,f,dx,f_tgl);
end

function adj(ℳ::Lorenz63Model,t,x,dx::AbstractVector)
    f_adj(t,x,dx)  = lorenz63_dxdt_adj(ℳ.σ,ℳ.β,ℳ.ρ,x,dx)
    f(t,x) = lorenz63_dxdt(ℳ.σ,ℳ.β,ℳ.ρ,x)

    return rungekutta2_adj(t,x,ℳ.dt,f,dx,f_adj);
end

function lorenz63_dxdt(σ,β,ρ,x)
    dxdt = similar(x)
    dxdt[1] = σ*(x[2]-x[1])
    dxdt[2] = x[1]*(ρ-x[3]) - x[2]
    dxdt[3] = x[1]*x[2] - β * x[3]
    return dxdt
end

function lorenz63_dxdt_tgl(σ,β,ρ,x,dx)
    Ddxdt = similar(x)
    Ddxdt[1] = σ*(dx[2]-dx[1])
    Ddxdt[2] = dx[1]*(ρ-x[3]) - x[1]*dx[3] - dx[2]
    Ddxdt[3] = dx[1]*x[2] + x[1]*dx[2] - β * dx[3]
    return Ddxdt
end


function lorenz63_dxdt_adj(σ,β,ρ,x,Ddxdt)
    dx = similar(x)
    #Ddxdt[1] = σ*(dx[2]-dx[1])
    dx[1] = -σ * Ddxdt[1]
    dx[2] = σ * Ddxdt[1]

    #Ddxdt[2] = dx[1]*(ρ-x[3]) - x[1]*dx[3] - dx[2]
    dx[1] += Ddxdt[2] * (ρ-x[3])
    dx[3] = -Ddxdt[2] * x[1]
    dx[2] += -Ddxdt[2]

    #Ddxdt[3] = dx[1]*x[2] + x[1]*dx[2] - β * dx[3]
    dx[1] += Ddxdt[3] * x[2]
    dx[2] += Ddxdt[3] * x[1]
    dx[3] += -Ddxdt[3] * β

    return dx
end


function rungekutta2(t,x,dt,f)
    k1 = dt * f(t,x);
    k2 = dt * f(t + dt/2,x + k1/2);

    xn = x + k2;
    return xn
end

function rungekutta2_tgl(t,x,dt,f,Dx,Df)
    k1 = dt * f(t,x);
    Dk1 = dt * Df(t,x,Dx);

    Dk2 = dt * Df(t + dt/2,x + k1/2,Dx + Dk1/2);
    Dxn = Dx + Dk2;
    return Dxn
end

function rungekutta2_adj(t,x,dt,f,Dxn,Df_adj)
    k1 = dt * f(t,x);

    Dx = Dxn;
    Dk2 = Dxn;

    Dtmp2 = Df_adj(t + dt/2, x+k1/2, dt*Dk2);
    Dx = Dx + Dtmp2;
    Dk1 = Dtmp2/2 ;
    Dx = Dx + dt*Df_adj(t,x,Dk1);

    return Dx
end


function rungekutta4(t,x,dt,f)
    k1 = dt * f(t,x);
    k2 = dt * f(t + dt/2,x + k1/2);
    k3 = dt * f(t + dt/2,x + k2/2);
    k4 = dt * f(t + dt  ,x + k3/2);
    xn = x + k1/6 + k2/3 + k3/3 + k4/6;
    return xn
end
