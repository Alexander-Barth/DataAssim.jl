
# nmax: total number of integration of the model
# x of size n x (nmax+1)

"""
    x,P = KalmanFilter(xi,Pi,ℳ,Q,yo,R,H,nmax,no)

Kalman Filter with the model `ℳ` and `nmax` time-steps starting at the
initial condition `xi` and error covariance `Pi`.
Observations `yo` (and error covariance `R`) at the time steps given in `no` are
assimilated with the observation operator `H`.
"""
function KalmanFilter(xi,Pi,ℳ,Q,yo::AbstractMatrix,R::AbstractMatrix,H::AbstractMatrix,nmax,no)

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


function KalmanFilter(xi,Pi,ℳ,Q,yo::Function,R::Function,𝓗::AbstractModel,nmax,no)
    T = eltype(xi)
    KalmanFilter(
        xi,Pi,ℳ,Q,
        VectorFun(T,length(no),yo),
        VectorFun(T,length(no),R),
        𝓗::AbstractModel,nmax,no)
end

function KalmanFilter(xi,Pi,ℳ,Q,yo::AbstractVector,R::AbstractVector,𝓗::AbstractModel,nmax,no)
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
            xn = @view x[:,n]
            Rn = R[obsindex]
            yon = yo[obsindex]
            HP = zeros(eltype(P),length(yon),length(xi))
            HPH = zeros(eltype(P),length(yon),length(yon))
            for i = 1:length(xi)
                HP[:,i] = tgl(𝓗,obsindex,xn,(@view P[:,i,n]))
            end
            for i = 1:length(yon)
                HPH[i,:] = tgl(𝓗,obsindex,xn,(@view HP[i,:]))
            end
            Hx = 𝓗(obsindex,xn)

            # assimilation
            K = HP' * inv(HPH + Rn)
            x[:,n]  = x[:,n] + K * (yo[obsindex] - Hx)
            P[:,:,n] = P[:,:,n] - K*HP

            obsindex = obsindex + 1
        end
    end

    return x,P
end
