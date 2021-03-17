function vout = softthr(vin, n, tau)

vin = reshape(vin, n, n);
[Uv, Sv, Vv] = svd(vin);
snew = max(diag(Sv) - tau,0);
vout = Uv*diag(snew)*Vv';
vout = vout(:);

end