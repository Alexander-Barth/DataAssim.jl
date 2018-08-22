# nmax: total number of integration of the model
# x of size n x (nmax+1)

function [x,P] = KalmanFilter(xi,Pi,model,Q,yo,R,H,nmax,no)

x = zeros(size(xi,1),size(yo,2));
P = zeros(size(xi,1),size(xi,1),size(yo,2));
obsindex = 1;

for n=1:nmax+1
    if n == 1
        # initialization
        x[:,1] = xi;
        P[:,:,1] = Pi;
    else
        # time integration
        x[:,n] = model.fun(n-1,x[:,n-1],0);
        P[:,:,n] = model.tgl(n-1,x[:,n-1],model.tgl(n-1,x[:,n-1],P[:,:,n-1])') + Q;
    end

    if obsindex <= length(no) && n == no(obsindex)
        # assimilation
        K = P[:,:,n]*H'*inv(H*P[:,:,n]*H' + R);
        x[:,n]  = x[:,n] + K * (yo[:,obsindex] - H*x[:,n]);
        P[:,:,n] = P[:,:,n] - K*H*P[:,:,n];

        obsindex = obsindex + 1;
    end
end
