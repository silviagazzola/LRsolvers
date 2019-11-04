function [Xall, RelRes, RelErr] = RS_GMRES_LRP(A,b,cycles,m,x0,truncbasis,truncsol,xex)

% input:
% iftrun - 1 for computing truncated SVD of V and X, 2 for doing singular 
%           value shrinkage (as in SVT)
% tau - shrinkage parameter (as in SVT)
% reg - regularization parameter for hybrid FGMRES; 'discrep' for discrepancy 
%           principle based parameter selection, or other numerical values
% eta, nl - used for discrepancy principle; eta is slightly larger than 1,
%           nl is noise level

N = length(x0);
totIt = cycles*m;
Xall = zeros(N,cycles);
RelErr = zeros(totIt, 1);
RelRes = zeros(totIt, 1);
totIt_temp = 0;

for i = 1:cycles
    [xapprox, RelRes_temp, RelErr_temp] = GMRES_LRP(A,b,m,x0,truncbasis,truncsol,xex);
    
    Xall(:,totIt_temp+1:totIt_temp+m) = xapprox;
    RelErr(totIt_temp+1:totIt_temp+m) = RelErr_temp;
    RelRes(totIt_temp+1:totIt_temp+m) = RelRes_temp;
    
    totIt_temp = totIt_temp + m;
    x0 = xapprox(:,end);
end

end

function [xapprox, RelRes, RelErr] = GMRES_LRP(A,b,m,x0,truncbasis,truncsol,xex)

N = size(x0,1);
n = sqrt(N);
xapprox = zeros(N,m);
RelRes = zeros(m,1);
RelErr = zeros(m,1);
nb = norm(b);
nx = norm(xex);

Vmatr = zeros(N,m+1);
W = zeros(N,m);
res = b - A*x0;
Vmatr(:,1)=res/norm(res);

for k = 1:m
    w = A*Vmatr(:,k);
    W(:,k) = w;
    alpha = (Vmatr(:,1:k)'*Vmatr(:,1:k))\(Vmatr(:,1:k)'*w);
    vtemp = w -Vmatr(:,1:k)*alpha;
    vtemp = truncfcn(vtemp, n, truncbasis);
    Vmatr(:,k+1)=vtemp/norm(vtemp);
    
    beta = (W(:,1:k)'*W(:,1:k))\(W(:,1:k)'*res);
    xtemp = x0 + Vmatr(:,1:k)*beta;
    xapprox(:,k) = truncfcn(xtemp, n, truncsol);
    RelErr(k) = norm(xapprox(:,k)-xex)/nx;
    RelRes(k) = norm(A*xapprox(:,k)-b)/nb;
end

end



