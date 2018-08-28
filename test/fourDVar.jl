# nmax: total number of integration of the model
# x of size n x (nmax+1)

function FreeRun!(ℳ,model_fun,xi,Q,H,nmax,no,x,Hx)
    obsindex = 1;

    for n=1:nmax+1
        if n == 1
            x[:,1] = xi;
        else
            x[:,n] = model_fun(n-1,x[:,n-1],0);
        end

        if obsindex <= length(no) && n == no[obsindex]
            # extract observations
            Hx[:,obsindex] = H*x[:,n];
            obsindex = obsindex +1;
        end
    end
end

function FreeRun(ℳ,model_fun,xi,Q,H,nmax,no)
    T = eltype(xi)
    x = zeros(T,size(xi,1),nmax+1);
    Hx = zeros(T,size(H,1),length(no))
    FreeRun!(ℳ,model_fun,xi,Q,H,nmax,no,x,Hx)
    return x,Hx
end


function cost2(xi,Pi,ℳ,model_fun,xa,yo,R,H,nmax,no,x,Hx)
    FreeRun!(ℳ,model_fun,xa,[],H,nmax,no,x,Hx);

    # cost function
    tmp = x[:,1] - xa;
    J = tmp' * (Pi \ tmp);
    for i = 1:size(yo,2)
        tmp = yo[:,i] - Hx[:,i];
        J = J + tmp' * (R \ tmp);
    end

    return J
end


function gradient(xi,dx0,x,Pi,model_tgl,model_adj,yo,R,H,nmax,no)

    dx = zeros(size(xi,1),nmax+1);
    dx[:,1] = dx0;
    obsindex = length(no);

    for n=1:nmax
        dx[:,n+1] = model_tgl(n,x[:,n],dx[:,n]);
    end

    lambda = zeros(size(xi,1),nmax+2);
    for n=nmax+1:-1:1
        lambda[:,n] = model_adj(n,x[:,n],lambda[:,n+1]);

        if obsindex > 0 && n == no[obsindex]
            lambda[:,n] = lambda[:,n] + H'*inv(R)*(yo[:,obsindex] - H*(dx[:,n]+x[:,n] ));
            obsindex = obsindex - 1;
        end
    end

    grad = inv(Pi)*(xi - (dx[:,1]+x[:,1])) + lambda[:,1];
    grad = -2 * grad;

    return grad,lambda
end


function fourDVar(
    xi::AbstractVector,Pi,ℳ,model_fun,model_tgl,model_adj,yo,R,H,nmax,no;
    innerloops = 10,
    outerloops = 2,
    tol = 1e-5)

#function fourDVar(xi::AbstractVector,Pi,ℳ,model_fun,model_tgl,model_adj,yo,R,H,nmax,no)
#    innerloops = 10
#    outerloops = 2
#    tol = 1e-5

    xa = float(xi)
    x = zeros(size(xi,1),nmax+1);
    Hx = zeros(size(yo,1),nmax+1);

    J = zeros(outerloops)
    b = zeros(size(xi))

    for i=1:outerloops

        # run with non-linear model
        #[x,Hx] = FreeRun(ℳ,model_fun,xa,[],H,nmax,no);
        #J[i] = cost(xi,Pi,x,yo,R,Hx);
        #J[i] = cost(xa,Pi,x,yo,R,Hx);
        J[i] = cost2(xi,Pi,ℳ,model_fun,xa,yo,R,H,nmax,no,x,Hx);
        #J[i],x,Hx = Jfun(xa);

        # dx increment relative to xi

        grad(dx) = gradient(xi,dx,x,Pi,model_tgl,model_adj,yo,R,H,nmax,no)[1];
        b .= grad(zeros(size(xi)));
        function fun(dx,fdx)
            fdx[:] = b - grad(dx)
        end

        dxa,cgsuccess,niter = DIVAnd.conjugategradient(fun,b,tol=tol/norm(b),maxit=innerloops,x0=zeros(size(xi)))

        @debug begin
            tmp = zeros(size(xi))
            fun(dxa,tmp)
            tmp = tmp - b

            @show tmp,sum(tmp.^2)
        end

        # add increment to dxa
        xa .= xa + dxa;
        #xa =  dxa;
    end

    return xa,J#,Jfun
end


