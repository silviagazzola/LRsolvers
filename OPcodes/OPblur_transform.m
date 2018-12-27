function vout = OPblur_transform(vin, A, W, transp_flag)

% OPblur_transform The forward computation and its adjoint modeling an 
% image blurring process, where sparsity transform may be incorporated
%
% w = OPblur_transform(vin, A, W, tflag)
%
% Performs the forward computation or returns the problem dimensions.
%
% Input: vin    - the vector to be operated upon
%          A    - the blurring matrix (can be a PSF matrix)
%          W    - sparsity transform
%        tflag  - string that determines the computation:
%                 'notransp', 'transp' or 'size'

% Julianne Chung, Virginia Tech
% Silvia Gazzola, University of Bath
% June, 2018.

% This file extends the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

[n, m] = size(vin);
if m == 1
    n = sqrt(n);
end
if strcmpi(tflag,'size')
    vout(1) = n;
    vout(2) = n;
elseif strcmp(transp_flag, 'notransp')
    if isscalar(W)
        if W == 1
            vout = reshape(A*reshape(vin,n,n),n*n,1);
        end
    else
        vout = reshape(W*(A*reshape(W'*reshape(vin,n,n), n,n)), n*n,1);
    end
else
    if isscalar(W)
        if W == 1
            vout = reshape(A'*reshape(vin,n,n),n*n,1);
        end
    else
        vout = reshape(W*(A'*reshape(W'*reshape(vin,n,n), n,n)), n*n,1);
    end
end