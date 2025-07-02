function TwinExperiment(ℳ,xit,Pi,Q,R,H,nmax,
                        no,method)


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
else
  error("unkown method: $(method)");
end


return xt,xfree,xa,yt,yo,diag

end

