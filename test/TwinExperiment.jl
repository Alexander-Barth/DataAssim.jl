function [xt,xfree,xa,yt,yo,diag] = TwinExperiment(model,xit,Pi,Q,R,H,nmax,...
                                                  no,method,varargin)

Nens = 100;

for i=1:2:length(varargin)
    if strcmp(varargin{i},'Nens')
        Nens = varargin{i+1};
    else
        error(['unkown argument ' varargin{i}]);
    end
end


n = length(xit);
m = size(R,1);
diag = struct();

% true run
[xt,yt] = FreeRun(model,xit,Q,H,nmax,no);

% add perturbations to IC   
xi = xit + chol(Pi)*randn(n,1);

% add perturbations to obs
yo = zeros(m,length(no));
for i=1:length(no)
  yo(:,i) = yt(:,i) + chol(R)*randn(m,1);
end

% free run
[xfree,yfree] = FreeRun(model,xi,Q,H,nmax,no);

% assimilation

if strcmp(method,'4DVar')
  [xia,J,Jfun] = fourDVar(xi,Pi,model,yo,R,H,nmax,no,'outerloops',10);
  [xa,ya] = FreeRun(model,xia,Q,H,nmax,no);  
  Pa = [];
  diag.J = J;
  diag.Jfun = Jfun;
elseif strcmp(method,'KF')
  [xa,Pa] = KalmanFilter(xi,Pi,model,Q,yo,R,H,nmax,no);
  diag.Pa = Pa;
elseif strcmp(method,'ETKF') || strcmp(method,'CLEnKF') || ...
      strmatch('sangoma-',method) == 1
  Ei = repmat(xi,[1 Nens]) + chol(Pi)*randn(n,Nens);
  [E] = ETKF_mod(Ei,model,Q,yo,R,H,nmax,no,method);
  xa = permute(mean(E,2),[1 3 2]);
  diag.E = E;
elseif strcmp(method,'KPF')
  Ei = repmat(xi,[1 Nens]) + chol(Pi)*randn(n,Nens);
  [E] = kpf(Ei,model,Q,yo,R,H,nmax,no,method);
  xa = permute(mean(E,2),[1 3 2]);
  diag.E = E;
else
  error(['unkown method:' method]);
end





%plot(xtrue);