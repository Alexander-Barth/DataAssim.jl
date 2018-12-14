
# nmax: total number of integration of the model
# x of size n x (nmax+1)

"""
    x,P = KalmanFilter(xi,Pi,ℳ,Q,yo,R,H,nmax,no)

Kalman Filter with the model `ℳ` and `nmax` time-steps starting at the
initial condition `xi` and error covariance `Pi`.
Observations `yo` (and error covariance `R`) at the time steps given in `no` are
assimilated with the observation operator `H`.
"""
function KalmanFilter(xi,Pi,ℳ,Q,yo,R,H,nmax,no)

    x = zeros(size(xi,1),nmax+1);
    P = zeros(size(xi,1),size(xi,1),nmax+1);
    obsindex = 1;

    for n=1:nmax+1
        if n == 1
            # initialization
            x[:,1] = xi;
            P[:,:,1] = Pi;
        else
            # time integration
            x[:,n] = ℳ(n-1,x[:,n-1]);
            P[:,:,n] = tgl(ℳ,n-1,x[:,n-1],tgl(ℳ,n-1,x[:,n-1],P[:,:,n-1])') + Q;
        end

        if obsindex <= length(no) && n == no[obsindex]
            # assimilation
            K = P[:,:,n]*H'*inv(H*P[:,:,n]*H' + R);
            x[:,n]  = x[:,n] + K * (yo[:,obsindex] - H*x[:,n]);
            P[:,:,n] = P[:,:,n] - K*H*P[:,:,n];

            obsindex = obsindex + 1;
        end
    end

    return x,P
end
