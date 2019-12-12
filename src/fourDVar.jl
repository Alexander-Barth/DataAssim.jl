# nmax: total number of integration of the model
# x of size n x (nmax+1)

function FreeRun!(ℳ,xi,Q,H,nmax,no,x,Hx)
    obsindex = 1;

    for n=1:nmax+1
        if n == 1
            x[:,1] = xi;
        else
            x[:,n] = ℳ(n-1,x[:,n-1]);
        end

        if obsindex <= length(no) && n == no[obsindex]
            # extract observations
            Hx[:,obsindex] = H*x[:,n];
            obsindex = obsindex +1;
        end
    end
end

"""
    x,Hx = FreeRun(ℳ,xi,Q,H,nmax,no)

Performs a free-run with the model `ℳ` and `nmax` time-steps starting at the
initial condition `xi`. Observations at the time steps given in `no` are 
extracted with the observation operator `H`.
"""
function FreeRun(ℳ,xi,Q,H,nmax,no)
    T = eltype(xi)
    x = zeros(T,size(xi,1),nmax+1);
    Hx = zeros(T,size(H,1),length(no))
    FreeRun!(ℳ,xi,Q,H,nmax,no,x,Hx)
    return x,Hx
end


function costfun(xi,Pi,ℳ,xa,yo,R,H,nmax,no,x,Hx)
    FreeRun!(ℳ,xa,[],H,nmax,no,x,Hx);

    # cost function
    tmp = x[:,1] - xa;
    J = tmp' * (Pi \ tmp);
    for i = 1:size(yo,2)
        tmp = yo[:,i] - Hx[:,i];
        J = J + tmp' * (R \ tmp);
    end

    return J
end


function gradient(xi,dx0,x,Pi,ℳ,yo,R,H,nmax,no)

    dx = zeros(size(xi,1),nmax+1);
    dx[:,1] = dx0;
    obsindex = length(no);

    for n=1:nmax
        dx[:,n+1] = tgl(ℳ,n,x[:,n],dx[:,n]);
    end

    lambda = zeros(size(xi,1),nmax+2);
    for n=nmax+1:-1:1
        lambda[:,n] = adj(ℳ,n,x[:,n],lambda[:,n+1]);

        if obsindex > 0 && n == no[obsindex]
            lambda[:,n] = lambda[:,n] + H'*inv(R)*(yo[:,obsindex] - H*(dx[:,n]+x[:,n] ));
            obsindex = obsindex - 1;
        end
    end

    grad = inv(Pi)*(xi - (dx[:,1]+x[:,1])) + lambda[:,1];
    grad = -2 * grad;

    return grad,lambda
end

"""
    x,J = fourDVar(
            xi,Pi,ℳ,yo,R,H,nmax,no;
            innerloops = 10,
            outerloops = 2,
            tol = 1e-5)

Incremental 4D-Var with the model `ℳ` and `nmax` time-steps starting at the
initial condition `xi` and error covariance `Pi` with the specified numbers of inner 
and outer loops.
Observations `yo` (and error covariance `R`) at the time steps given in `no` are
assimilated with the observation operator `H`.
"""
function fourDVar(
    xi::AbstractVector,Pi,ℳ,yo,R,H,nmax,no;
    innerloops = 10,
    outerloops = 2,
    tol = 1e-5)

#function fourDVar(xi::AbstractVector,Pi,ℳ,yo,R,H,nmax,no)
#    innerloops = 10
#    outerloops = 2
#    tol = 1e-5

    xa = float(xi)
    x = zeros(size(xi,1),nmax+1);
    Hx = zeros(size(yo,1),nmax+1);

    J = zeros(outerloops)
    b = zeros(size(xi))

    for i=1:outerloops

        J[i] = costfun(xi,Pi,ℳ,xa,yo,R,H,nmax,no,x,Hx);

        # dx increment relative to xi

        grad(dx) = gradient(xi,dx,x,Pi,ℳ,yo,R,H,nmax,no)[1];
        b .= grad(zeros(size(xi)));

        #=
        function fun(dx,fdx)
            #@show "D"
            fdx[:] = b - grad(dx)
        end

        dxa,cgsuccess,niter = DIVAnd.conjugategradient(fun,b,tol=tol/norm(b),maxit=innerloops,x0=zeros(size(xi)))
        @debug begin
            tmp = zeros(size(xi))
            fun(dxa,tmp)
            tmp = tmp - b

            @show "DIVAnd.conjugategradient",sum(tmp.^2),niter
        end
        =#

        function fg!(F,G,x)
            @debug "call fg! $(F==nothing) $(G==nothing)"
            GG = grad(x)
            if G != nothing
                # code to compute gradient here
                # writing the result to the vector G
                G .= GG
            end
            if F != nothing
                value = 0.5 * (x ⋅ GG) + 0.5 * b ⋅ x
                return value
            end
        end

        result = Optim.optimize(Optim.only_fg!(fg!), zeros(size(b)), ConjugateGradient(),
                                Optim.Options(g_tol = tol,
                                              iterations = innerloops,
                                              allow_f_increases=true,
                                              store_trace = true,
                                              show_trace = false))
        dxa = result.minimizer

        @debug begin
            @show summary(result)
            tmp = zeros(size(xi))
            fg!(nothing,tmp,dxa)
            @show tmp

            @show innerloops
            @show Optim.g_converged(result)
            @show Optim.f_converged(result)
            @show "optim",sum(tmp.^2),Optim.iterations(result)
        end

        # add increment to dxa
        xa .= xa + dxa;
    end

    return xa,J#,Jfun
end


