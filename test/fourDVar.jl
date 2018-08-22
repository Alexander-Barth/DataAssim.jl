# nmax: total number of integration of the model
# x of size n x (nmax+1)

function FreeRun(model_fun,xi,Q,H,nmax,no)

    x = zeros(size(xi,1),nmax+1);
    obsindex = 1;
    Hx = zeros(size(H,1),length(no))

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

    return x,Hx
end


function cost2(xi,Pi,model_fun,xa,yo,R,H,nmax,no)
    x,Hx = FreeRun(model_fun,xa,[],H,nmax,no);

    # cost function
    tmp = x[:,1] - xa;
    J = tmp' * (Pi \ tmp);
    for i = 1:size(yo,2)
        tmp = yo[:,i] - Hx[:,i];
        J = J + tmp' * (R \ tmp);
    end

    return J,x,Hx
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

    #-2*(inv(Pi)*(xi - x0)  + H'*inv(R)*(yo - H*xi))
    return grad,lambda
end


function fourDVar(xi::AbstractVector,Pi,model_fun,model_tgl,model_adj,yo,R,H,nmax,no; innerloops = 10,
    outerloops = 2,
    tol = 1e-5)

    xa = xi
    x = zeros(size(xi,1),nmax+1);

    Jfun(xa) = cost2(xi,Pi,model_fun,xa,yo,R,H,nmax,no);
    J = zeros(outerloops)

    for i=1:outerloops

        # run with non-linear model
        #[x,Hx] = FreeRun(model_fun,xa,[],H,nmax,no);
        #J[i] = cost(xi,Pi,x,yo,R,Hx);
        #J[i] = cost(xa,Pi,x,yo,R,Hx);
        #[J[i],x,Hx] = cost2(xi,Pi,model,xa,yo,R,H,nmax,no);
        J[i],x,Hx = Jfun(xa);

        # dx increment relative to xi

        grad(dx) = gradient(xi,dx,x,Pi,model_tgl,model_adj,yo,R,H,nmax,no)[1];
        b = grad(zeros(size(xi)));
        function fun(dx,fdx)
            fdx[:] = b - grad(dx);
        end

        dxa,cgsuccess,niter = DIVAnd.conjugategradient(fun,b,tol=tol/norm(b),maxit=innerloops,x0=zeros(size(xi)))

        @debug begin
            tmp = zeros(size(xi))
            fun(dxa,tmp)
            tmp = tmp - b

            @show tmp,sum(tmp.^2)
        end

        # add increment to dxa
        xa = xa + dxa;
        #xa =  dxa;
    end

    return xa,J,Jfun
end


