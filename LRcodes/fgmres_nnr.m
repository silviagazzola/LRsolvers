function [xapprox,RelErr,RelRes,lambdahis,Zmatr,Vmatr,H] = fgmres_nnr(A, b, xex, x0,parameters)

m = parameters.maxIt;
reg = parameters.reg;
regmat = parameters.regmat;
svdbasis = parameters.svdbasis;
eta = parameters.eta;
nl = parameters.nl;
p = parameters.p;

N = length(b); n = sqrt(N);
nim = sqrt(N);

xapprox = zeros(N,m);
RelRes = zeros(m,1);
RelErr = zeros(m,1);

nb = norm(b);
nx = norm(xex);


B = reshape(b,nim,nim);
[UB,~,VB] = svd(B);

UX = UB; VX = VB; % in the case x0 = 0;
SX = eye(nim); 
%precX = SX(:);
sx = diag(SX); 
precX = kron(ones(n,1),sx);
precX = precX.^((2-p)/2);

Vmatr = zeros(N,m+1);
Zmatr = zeros(N,m);
res = b - A*x0;
Vmatr(:,1)=res/norm(res);
nr = norm(res);

etaeps = eta*nl;

if strcmp(reg,'discrep')
    lambda = 1e8;
    % lambda = 1;
elseif isnumeric(reg)
    lambda = reg;
end
lambdahis = ones(m,1);

for k = 1:m
    
%%%%
  vtemp = Vmatr(:,k);
  vtemp = reshape(vtemp, nim, nim);
  [UV, ~, VV] = svd(vtemp);
  if svdbasis == 1
      vtemp = UV'*(vtemp*VV);
      vtemp = precX.*vtemp(:);
      vtemp = reshape(vtemp, nim, nim); 
      vtemp = UV*(vtemp*VV');
  elseif svdbasis == 2
      vtemp = UX'*(vtemp*VX);
      vtemp = precX.*vtemp(:);
      vtemp = reshape(vtemp, nim, nim); 
      vtemp = UV*(vtemp*VV');
  elseif svdbasis == 3
      vtemp = UB'*(vtemp*VB);
      vtemp = precX.*vtemp(:);
      vtemp = reshape(vtemp, nim, nim); 
      vtemp = UV*(vtemp*VV');
  elseif svdbasis == 4
      vtemp = UV'*(vtemp*VV);
      vtemp = precX.*vtemp(:);
      vtemp = reshape(vtemp, nim, nim); 
      vtemp = UX*(vtemp*VX');
  elseif svdbasis == 5
      vtemp = UX'*(vtemp*VX);
      vtemp = precX.*vtemp(:);
      vtemp = reshape(vtemp, nim, nim); 
      vtemp = UX*(vtemp*VX');
  elseif svdbasis == 6
      vtemp = UB'*(vtemp*VB);
      vtemp = precX.*vtemp(:);
      vtemp = reshape(vtemp, nim, nim); 
      vtemp = UX*(vtemp*VX');
%   elseif svdbasis == 7
%       vtemp = UV'*(vtemp*VV);
%       vtemp = precX.*vtemp(:);
%       vtemp = reshape(vtemp, nim, nim); 
%       vtemp = UB*(vtemp*VB');
%   elseif svdbasis == 8
%       vtemp = UX'*(vtemp*VX);
%       vtemp = precX.*vtemp(:);
%       vtemp = reshape(vtemp, nim, nim); 
%       vtemp = UB*(vtemp*VB');
%   elseif svdbasis == 9
%       vtemp = UB'*(vtemp*VB);
%       vtemp = precX.*vtemp(:);
%       vtemp = reshape(vtemp, nim, nim); 
%       vtemp = UB*(vtemp*VB');
  end
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

        if strcmp(regmat,'I')
            c = UH'*d;

            if k == 1
                SH = SH(1,1);
            else
                SH = diag(SH);
            end

            Filt = SH.^2 + lambda;
            yj = VH*((SH.*c(1:k))./Filt);
            
        elseif strcmp(regmat,'R')
            [~,R] = qr(Zmatr(:,1:k),0);
            yj = [H(1:k+1,1:k);sqrt(lambda)*R]\[d;zeros(k,1)];
        end
        
        RelRes(k) = norm(H(1:k+1,1:k)*yj - d)/nb;
        
        if strcmp(reg,'discrep')
            beta = 1/lambda;
            SH = [diag(SH); zeros(1,k)];
            matrtemp = (beta*(SH*SH') + eye(k+1));
            zetabeta = matrtemp\c;
            wbeta = matrtemp\zetabeta; 
            
            f = (RelRes(k))^2;
            f1 = 2/beta*zetabeta'*(wbeta - zetabeta); f1 = f1/(nr^2);
            
            % f1 = ...
            beta=beta-(f-nl^2)/(f1); 
            lambda = 1/beta;
            
%             lambda = abs((etaeps-gmres_res)/(RelRes(k)-gmres_res))*lambda;
            lambdahis(k) = lambda; 
      
        end
        
    end
    xj = x0 + Zmatr(:,1:k)*yj;
    
    xapprox(:,k) = xj(:);
    RelErr(k) = norm(xj-xex)/nx;
    
    Xtemp = reshape(xj,nim,nim);
    [UX,SX,VX] = svd(Xtemp);
%     precX = SX(:);
%     precX = precX.^((2-p)/2);
    sx = diag(SX); 
    precX = kron(ones(n,1),sx);
    precX = precX.^((2-p)/2);

end

RelRes = [nr;RelRes];
RelErr = [norm(xex-x0)/nx;RelErr];