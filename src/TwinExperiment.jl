function TwinExperiment(ℳ,xit,Pi,Q,R,H,nmax,
                        no,method; Nens = 100)


n = length(xit);
diag = Dict();

# true run
xt,yt = FreeRun(ℳ,xit,Q,H,nmax,no);

# add perturbations to IC
xi = xit + cholesky(Pi).U * randn(n)

T = eltype(xit)
# add perturbations to obs
yo = Vector{Vector{T}}(undef,length(no))
for i = 1:length(no)
  m = length(yt[i])
  yo[i] = yt[i] + cholesky(R[i]).U * randn(m);
end

# free run
xfree,yfree = FreeRun(ℳ,xi,Q,H,nmax,no);

# assimilation

if method == "4DVar"
  xia,J = fourDVar(xi,Pi,ℳ,yo,R,H,nmax,no,outerloops = 10);
  xa,ya = FreeRun(ℳ,xia,Q,H,nmax,no);
  Pa = [];
  diag[:J] = J;
  #diag[:Jfun] = Jfun;
elseif method == "KF"
  xa,Pa = KalmanFilter(xi,Pi,ℳ,Q,yo,R,H,nmax,no);
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

