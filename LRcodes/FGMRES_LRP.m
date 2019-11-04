function [xapprox, RelRes, RelErr, lambdahis,Zmatr, Vmatr, H] = FGMRES_LRP(A,b,m,x0,truncbasis,truncsol,xex,tau,iftrun,reg,eta,nl)

% input:
% iftrun - 1 for computing truncated SVD of V and X, 2 for doing singular 
%           value shrinkage (as in SVT)
% tau - shrinkage parameter (as in SVT)
% reg - regularization parameter for hybrid FGMRES; 'discrep' for discrepancy 
%           principle based parameter selection, or other numerical values
% eta, nl - used for discrepancy principle; eta is slightly larger than 1,
%           nl is noise level



N = size(x0,1);
n = sqrt(N);
xapprox = zeros(N,m);
RelRes = zeros(m,1);
RelErr = zeros(m,1);
nb = norm(b);
nx = norm(xex);

Vmatr = zeros(N,m+1);
Zmatr = zeros(N,m);
res = b - A*x0;
Vmatr(:,1)=res/norm(res);
nr = norm(res);

etaeps = eta*nl;

if strcmp(reg,'discrep')
    lambda = 1;
elseif isnumeric(reg)
    lambda = reg;
end
lambdahis = ones(m,1);


for k = 1:m
  
    z = Vmatr(:,k);

    if iftrun == 1
        z = truncfcn(z, n, truncbasis);
    elseif iftrun == 2
        z = softthr(z, n,tau);
    end
    
    Zmatr(:,k) = z;  
    
    w = A*Zmatr(:,k);
    for i=1:k
        H(i,k)=w'*Vmatr(:,i);
        w=w-H(i,k)*Vmatr(:,i);
    end
    H(k+1,k)=norm(w);
    if H(k+1,k)<1e-10
        disp('Breakdown of GMRES')
        return
    else
        Vmatr(:,k+1)=w/H(k+1,k);
        d=[nr;zeros(k,1)];
        
        [UH,SH,VH] = svd(H(1:k+1,1:k));
        gmres_res = abs((UH(:,k+1)'*d))/nb;
        
        c = UH'*d;
        
        if k == 1
            SH = SH(1,1);
        else
            SH = diag(SH);
        end
        
        Filt = SH.^2 + lambda;
        yj = VH*((SH.*c(1:k))./Filt);
        
        RelRes(k) = norm(H(1:k+1,1:k)*yj - d)/nb;
        
        if strcmp(reg,'discrep')

            lambda = abs((etaeps-gmres_res)/(RelRes(k)-gmres_res))*lambda;
            lambdahis(k) = lambda; 
         
        end
    end
    xj = x0 + Zmatr(:,1:k)*yj;
    
    if iftrun == 1
        xj = truncfcn(xj, n, truncsol);
    elseif iftrun == 2
        xj = softthr(xj,n,tau);
    end
    
    xapprox(:,k) = xj(:);
    RelErr(k) = norm(xj-xex)/nx;
end



