using DataAssim
using Base.Test

# write your own tests here



# number of elements in the state vector
n = 10;
# ensemble size
N = 3;
# number of observations
m = 5;

# if debug is one, then internal checks are activated
debug = false;
debug = true;
tol = 1e-10;

y = randn(m,1);
Xf = randn(n,N);

H = randn(m,n);

y = 1:m
Xf = reshape(sin(3*(1:(n*N))),n,N);
H = reshape(1:(m*n),m,n);

R = 2*eye(m,m);

xf = mean(Xf,2);
Xfp = Xf - repmat(xf,1,N);


Pf = (Xfp * Xfp') / (N-1);
K = Pf * H'*inv(H*Pf*H' + R);
Pa_check = Pf - K*H*Pf;
xa_check = xf + K*(y - H*xf);

method = ["EnSRF","EAKF","ETKF","ETKF2","SEIK","ESTKF","serialEnSRF"]
#method = ["EnSRF","serialEnSRF","ETKF","ETKF2","EAKF","SEIK","ESTKF"]
method = ["EnSRF","EAKF","ETKF","ETKF2","SEIK","ESTKF"]

print("Non-serial algorithms (all observations at once) with H\n")

for i = 1:length(method)
  Xa,xa = sangoma_ensemble_analysis(Xf,H,y,R,method[i],
                                    debug=debug, tolerance=tol);
  Xap = Xa - repmat(xa,1,N);

  # check results
  sangoma_check(xa,xa_check,method[i] * "-analysis",tol)

  sangoma_check(mean(Xa,2),xa_check,
                method[i] * "-analysis ensemble mean",tol);

  sangoma_check((Xap * Xap') / (N-1),Pa_check,
                method[i] * "-analysis ensemble variance",tol)
end

# #
method = ["EnSRF","EAKF","ETKF","ETKF2","SEIK","ESTKF"]
#method = ["EnSRF"]


print("Non-serial algorithms (all observations at once) with HXf\n")

for i = 1:length(method)
  Xa,xa = sangoma_ensemble_analysis(Xf,[],y,R,method[i],
                                    debug=debug, tolerance=tol,HXf=H*Xf)
  Xap = Xa - repmat(xa,1,N)

  # check results
  sangoma_check(xa,xa_check,method[i] * "-analysis",tol)

  sangoma_check(mean(Xa,2),xa_check,
                method[i] * "-analysis ensemble mean",tol);

  sangoma_check((Xap * Xap') / (N-1),Pa_check,
                method[i] * "-analysis ensemble variance",tol)
end


#   method = {"EnSRF","ETKF","SEIK","ESTKF","serialEnSRF"};
#   print("Serial algorithms (one observation after the other)")


# for i = 1:length(method)
#   [Xa,xa] = sangoma_ensemble_analysis(Xf,H,y,R,method[i],
#       "debug",0,"tolerance",tol,"serial",1);
#   Xap = Xa - repmat(xa,[1 N]);

#   % check results
#   sangoma_check(xa,xa_check,method[i] * "-analysis",tol)

#   sangoma_check(mean(Xa,2),xa_check,
#                 method[i] * "-analysis ensemble mean",tol);

#   sangoma_check((Xap * Xap') / (N-1),Pa_check,
#                 method[i] * "-analysis ensemble variance",tol)
# end

