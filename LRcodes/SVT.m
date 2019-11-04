function [Xapprox,relerr,relres] = SVT(A,b,xtrue,ifSVT,m,tau,delta)

n1 = length(b);
n2 = length(xtrue);
n = sqrt(n2);

Xapprox = zeros(n2,m);
relerr = zeros(m,1);
relres = relerr;

y = zeros(n1,1);

for i = 1:m
    
        ATy = A'*y;
        
        
        if ifSVT == 1
            [U,S,V] = svd(reshape(ATy,n,n));
            snew = max(diag(S) - tau,0);
            X = U*diag(snew)*V';
            x = X(:);
        elseif ifSVT == 0
            x = ATy;
        else 
            [U,S,V] = svd(reshape(ATy,n,n));
            snew = diag(S);
            snew(snew<tau) = 0;
            X = U*diag(snew)*V';
            x = X(:);
        end
        
        
        r = b - A*x;
        y = y + delta*r;
        
    
    Xapprox(:,i) = x;
    relerr(i) = norm(x-xtrue)/norm(xtrue);
    relres(i) = norm(r)/norm(b);
end

