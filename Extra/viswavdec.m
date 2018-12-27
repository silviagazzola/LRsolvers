function im = viswavdec(C, A, p)
%
% This function takes the wavelet coefficients out of C (that were created
% using A) and presents a visually pleasing image.
% p is a scalar for better visualization (p=2 is default)

if nargin <3
  p = 2;
end

C = C(:);
wlev = A.wlevels;

[~,S] = wavedec2(rand(A.imsize),wlev,A.wname);

[cA] = appcoef2(C,S,A.wname,wlev); % approximation coefficients
im = cA/p;
for i = wlev:-1:1
  [cH,cV,cD] = detcoef2('all',C,S,i); 
  im = [im, cH;cV,cD]/p;
end
im = (p^wlev)*abs(im);
