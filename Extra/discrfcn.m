function discrval = discrfcn(l, A, RegM, b, nnoise)

n = size(A,2);
if n == 1
    discrval = 0;
else
    xl = [A; l*RegM]\[b; zeros(size(RegM,1),1)];
    discrval = (norm(A*xl -b)/norm(b))^2 - nnoise^2;
end