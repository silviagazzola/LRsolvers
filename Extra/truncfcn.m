function vout = truncfcn(vin, n, truncind)

vin = reshape(vin, n, n);
[Uv, Sv, Vv] = svd(vin);
Uv = Uv(:,1:truncind);
Sv = Sv(1:truncind,1:truncind);
Vv = Vv(:,1:truncind);
vout = Uv*Sv*Vv';
vout = vout(:);