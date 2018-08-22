function TwinExperiment(model_fun,model_tgl,model_adj,xit,Pi,Q,R,H,nmax,
                        no,method; Nens = 100)


n = length(xit);
m = size(R,1);
diag = Dict();

# true run
xt,yt = FreeRun(model_fun,xit,Q,H,nmax,no);

# add perturbations to IC
xi = xit + cholesky(Pi).U * randn(n)

# add perturbations to obs
yo = zeros(m,length(no));
for i = 1:length(no)
  yo[:,i] = yt[:,i] + cholesky(R).U * randn(m,1);
end

# free run
xfree,yfree = FreeRun(model_fun,xi,Q,H,nmax,no);

# assimilation

if method == "4DVar"
  xia,J,Jfun = fourDVar(xi,Pi,model_fun,model_tgl,model_adj,yo,R,H,nmax,no,outerloops = 10);
  xa,ya = FreeRun(model_fun,xia,Q,H,nmax,no);
  Pa = [];
  diag[:J] = J;
  diag[:Jfun] = Jfun;
elseif method == "KF"
  xa,Pa = KalmanFilter(xi,Pi,model_fun,model_tgl,Q,yo,R,H,nmax,no);
  diag[:Pa] = Pa;
# elseif strcmp(method,"ETKF") || strcmp(method,"CLEnKF") || ...
#       strmatch("sangoma-",method) == 1
#   Ei = repmat(xi,[1 Nens]) + cholesky(Pi).U * randn(n,Nens);
#   [E] = ETKF_mod(Ei,model,Q,yo,R,H,nmax,no,method);
#   xa = permute(mean(E,2),[1 3 2]);
#   diag.E = E;
# elseif strcmp(method,"KPF")
#   Ei = repmat(xi,[1 Nens]) + cholesky(Pi).U * randn(n,Nens);
#   [E] = kpf(Ei,model,Q,yo,R,H,nmax,no,method);
#   xa = permute(mean(E,2),[1 3 2]);
#   diag.E = E;
else
  error("unkown method: $(method)");
end


return xt,xfree,xa,yt,yo,diag

end

