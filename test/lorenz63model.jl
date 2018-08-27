mutable struct Lorenz63Model{T} <: AbstractModel
    dt::T
    sigma::T
    beta::T
    rho::T
end


Lorenz63Model(dt,sigma=10.,rho = 28,beta = 8/3) = Lorenz63Model(dt,rho,beta)

function (M::Lorenz63Model)(t,x,eta = 0)
    return rungekutta2(t,x,M.dt,(t,x) -> lorenz63_dxdt(M.sigma,M.beta,M.rho,x)) + eta
end

# function tgl(M::Lorenz63Model,t,dx)
#     M.M*dx
# end

# function adj(M::Lorenz63Model,t,dx)
#     M.M'*dx
# end



function lorenz63_dxdt(sigma,beta,rho,x)
    dxdt = similar(x)
    dxdt[1] = sigma*(x(2)-x(1))
    dxdt[1] = x(1)*(rho-x(3)) - x(2)
    dxdt[3] = x(1)*x(2) - beta * x(3)
    return dxdt
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


# function Dx = rungekutta2_adj_orig(t,x,dt,f,Dxn,Df_adj)

# k1 = dt * f(t,x);
# k2 = dt * f(t + dt/2,x + k1/2);
# xn = x + k2;

# # Dk1 = dt * Df(t,x,Dx);
# # # Dk1  = [0   dt*Df] Dk1
# # # Dx   = [0   1    ] Dx



# # Dtmp2 = Dx + Dk1/2;
# # # Dtmp2    [0  1 1/2] Dtmp2
# # # Dx    =  [0  1  0 ] Dx
# # # Dk1      [0  0   1] Dk1

# # Dtmp = Df(t + dt/2,x + k1/2,Dtmp2)
# # # Dtmp  = [0   Df] Dtmp
# # # Dtmp2   [0   1 ] Dtmp2


# # Dk2 = dt * Dtmp;
# # # Dk2    = [0  dt] Dk2
# # # Dtmp     [0  1 ] Dtmp

# # Dxn = Dx + Dk2;
# # #[0 1 1; 
# # # 0 1 0
# # # 0 0 1]

# Dx = 0;
# Dk1 = 0;
# Dk2 = 0;
# Dtmp = 0;
# Dtmp2 = 0;

# Dx = Dx + Dxn;
# Dk2 = Dk2 + Dxn;
# Dxn = 0;

# # Dk2    = [0   0] Dk2
# # Dtmp     [dt  1 ] Dtmp

# Dtmp = dt*Dk2 + Dtmp;
# Dk2 = 0;

# # Dtmp  = [0     0 ] Dtmp
# # Dtmp2   [Df'   1 ] Dtmp2

# Dtmp2 = Df_adj(t + dt/2, x+k1/2, Dtmp) + Dtmp2;

# # Dtmp2    [0  0  0 ] Dtmp2
# # Dx    =  [1  1  0 ] Dx
# # Dk1      [1/2 0 1 ] Dk1

# Dx = Dx + Dtmp2;
# Dk1 = Dk1 + Dtmp2/2 ;
# Dtmp2 = 0;

# # Dk1  = [0        0 ] Dk1
# # Dx   = [dt*Df'   1 ] Dx

# Dx = Dx + dt*Df_adj(t,x,Dk1);
# Dk1 = 0;

# end




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
