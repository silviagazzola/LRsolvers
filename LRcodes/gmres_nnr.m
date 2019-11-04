function [xapprox,RelErr,RelRes,lambdahis,Zmatr,Vmatr,H] = gmres_nnr(A, b, xex, x0,parameters)

m = parameters.maxIt;
reg = parameters.reg;
eta = parameters.eta;
nl = parameters.nl;
p = parameters.p;
UX = parameters.svd.U;
VX = parameters.svd.V;
SX = parameters.svd.S;
% precX = SX(:);
% precX = precX.^((2-p)/2);

N = length(x0); n = sqrt(N);
nim = sqrt(N);

sx = diag(SX); 
precX = kron(ones(n,1),sx);
precX = precX.^((2-p)/2);

xapprox = zeros(N,m);
RelRes = zeros(m,1);
RelErr = zeros(m,1);

nb = norm(b);
nx = norm(xex);


% B = reshape(b,nim,nim);
% [UB,~,VB] = svd(B);
% 
% UX = UB; VX = VB; % in the case x0 = 0;
% SX = eye(nim); precX = SX(:);

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
elseif isstruct(reg)
    % rule = reg.rule;
    lambda = reg.value;
    stagn = 0;
end
lambdahis = ones(m+1,1); % history of the lambda's

for k = 1:m
    
%%%%
  vtemp = Vmatr(:,k);
  vtemp = reshape(vtemp, nim, nim);
  
  vtemp = UX'*(vtemp*VX);
  vtemp = precX.*vtemp(:);
  vtemp = reshape(vtemp, nim, nim); 
  vtemp = UX*(vtemp*VX');
  Zmatr(:,k) = vtemp(:);
  
  
%%%%

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
        
        if strcmp(reg,'discrep') || (isstruct(reg) && strcmp(reg.rule,'discrep'))
            lambda = abs((etaeps-gmres_res)/(RelRes(k)-gmres_res))*lambda;
            % lambdahis(k+1) = lambda; 
        elseif (isstruct(reg) && strcmp(reg.rule,'discrepvar')) % && k>1
%             RegParamk = fzero(@(l)discrfcn(l, Mk, ZRksq, rhsk, eta*NoiseLevel), [0, 1e10]);
%             RegParamVect(k) = RegParamk;
            betaold = 1/lambda;
            SH = [diag(SH); zeros(1,k)];
            matrtemp = (betaold*(SH*SH') + eye(k+1));
            zetabeta = matrtemp\c;
            wbeta = matrtemp\zetabeta; 
            
            f = (RelRes(k))^2;
            f1 = 2/betaold*zetabeta'*(wbeta - zetabeta); f1 = f1/(nr^2);
            
            % f1 = ...
            beta=betaold-(f-nl^2)/(f1); 
            lambda = 1/beta;
            if lambda <= 1e-8
                lambda = 1/betaold;
                if lambda <=1e-14
                    stagn = stagn + 1;
                end
            end
            if stagn > 4
                break
            end
        end
        lambdahis(k+1) = lambda; 
    end
        
xj = x0 + Zmatr(:,1:k)*yj;
    
xapprox(:,k) = xj(:);
RelErr(k) = norm(xj-xex)/nx;
    
    % Xtemp = reshape(xj,nim,nim);
    
end

