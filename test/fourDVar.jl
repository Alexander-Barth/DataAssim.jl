function cost2(xi,Pi,model,xa,yo,R,H,nmax,no)

    x,Hx = FreeRun(model,xa,[],H,nmax,no);

    # cost function
    tmp = x[:,1] - xa;
    J = tmp' * (Pi \ tmp);

    for i = 1:size(yo,2)
        tmp = yo[:,i] - Hx[:,i];
        J = J + tmp' * (R \ tmp);
    end

    return J,x,Hx
end


function gradient(xi,dx0,x,Pi,model,yo,R,H,nmax,no)

    dx = zeros(size(xi,1),nmax+1);
    dx[:,1] = dx0;
    obsindex = length(no);

    for n=1:nmax
        dx(:,n+1) = model.tgl(n,x(:,n),dx(:,n));
    end

    lambda = zeros(size(xi,1),nmax+2);
    for n=nmax+1:-1:1
        lambda(:,n) = model.adj(n,x(:,n),lambda(:,n+1));

        if obsindex > 0 && n == no(obsindex)
            lambda(:,n) = lambda(:,n) + H'*inv(R)*(yo(:,obsindex) - H*(dx(:,n)+x(:,n) ));
            obsindex = obsindex - 1;
        end
    end

    #grad = inv(Pi)*(xi - x[:,1]) + lambda[:,1];
    grad = inv(Pi)*(xi - (dx[:,1]+x(:,n))) + lambda[:,1];
    grad = -2 * grad;

    #-2*(inv(Pi)*(xi - x0)  + H'*inv(R)*(yo - H*xi))

    return grad,lambda
end


function fourDVar(xi,Pi,model,yo,R,H,nmax,no; innerloops = 10,
    outerloops = 2,
    tol = 1e-5)

    prop = varargin;
    for i=1:2:length(prop)
        if strcmp(prop{i},'innerloops')
            innerloops  = prop{i+1};
        elseif strcmp(prop{i},'outerloops')
            outerloops  = prop{i+1};
        elseif strcmp(prop{i},'tol')
            tol  = prop{i+1};
        else
        end
    end

    xa = xi;
    x = zeros(size(xi,1),nmax+1);

    Jfun = @(xa) cost2(xi,Pi,model,xa,yo,R,H,nmax,no);

    for i=1:outerloops

        # run with non-linear model
        #[x,Hx] = FreeRun(model,xa,[],H,nmax,no);
        #J(i) = cost(xi,Pi,x,yo,R,Hx);
        #J(i) = cost(xa,Pi,x,yo,R,Hx);
        #[J(i),x,Hx] = cost2(xi,Pi,model,xa,yo,R,H,nmax,no);
        [J(i),x,Hx] = Jfun(xa);

        # dx increment relative to xi

        grad = @(dx) gradient(xi,dx,x,Pi,model,yo,R,H,nmax,no);
        b = grad(zeros(size(xi)));
        fun = @(dx) b - grad(dx);

        [dxa] = conjugategradient(fun,b,'tol',tol,'maxit',innerloops,'x0',zeros(size(xi)));

        sum((fun(dxa)-b).^2)
        #  assert(sum((fun(dxa)-b).^2) < 10*tol^2)

        # add increment to dxa
        xa = xa + dxa;
        #xa =  dxa;
    end

    return xa,J,Jfun
end


