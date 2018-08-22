# nmax: total number of integration of the model
# x of size n x (nmax+1) 

function [x,Hx] = FreeRun(model,xi,Q,H,nmax,no)

x = zeros(size(xi,1),nmax+1);
obsindex = 1;

for n=1:nmax+1
    if n == 1
        x(:,1) = xi;
    else
        x[:,n] = model.fun(n-1,x(:,n-1),0);
    end
    
    if obsindex <= length(no) && n == no(obsindex)
      # extract observations
      Hx(:,obsindex) = H*x[:,n];        
      obsindex = obsindex +1;
    end
end
