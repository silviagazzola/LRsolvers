function E = TikOptParam(lambda, H, L, d, x0, Z, W, xex)
%
% Function used to compute the optimal regularization parameter for the
% projected problem
%
% Input:
%   lambda - regularization parameter
%   H - matrix of projected problem
%   L - regularization matrix for projected problem
%   d - right hand side for projected problem
%   x0 - initital vector to update
%   Z - basis vectors
%   W - transfomration matrix
%   xex - exact solution
%
% Output:
%   E = ||x_j - x_true||/||x_true||
%
% J. Chung, 2/8/18


p = 2; %2-norm of error
j = size(L,1);

% solve projected problem
if lambda<=eps
  yj = H\d;
else
  Mj = [H; lambda*L];        
  % Mj = [H; lambda*L];        
  cj = [d; zeros(j,1)]; 
  yj = Mj \ cj;
%   Mj = H'*H + lambda*L'*L;        
%   cj = H'*d;
%   yj = Mj \ cj;
end

% Back into original space
xj = x0(:) + Z*yj;

% Undo transformation
xj = W'*xj; 

% Relative error
% E = norm(xj(:) - xex(:))/norm(xex(:));
E = norm(xj(:) - xex(:),p)/norm(xex(:),p);
% E = norm(yj(:) - Z \ (W*xex(:)-x0(:)),p);