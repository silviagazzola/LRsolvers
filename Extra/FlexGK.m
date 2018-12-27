function [U, M, T, V, z] = FlexGK(A, U, M, T, V, L)
%
%     [U, M, T, V, z] = FlexGK(A, U, M, T, V, L)
%
%  Perform one step of (inexact/flexible) Golub Kahan with
%  reorthogonalization.
%
% Input:
%          A - matrix A and A-transpose functions
%       U, V - accumulation of orthonormal vectors
%          M - k x k upper triangular matrix
%          T - (k+1) x k upper Hessenberg matrix
%          L - Changing preconditioner matrix
%
% Output:
%       U, V - updated matrix or orthonormal columns
%       M, T - updated matrices
%          z - new z vector
%

% Julianne Chung, Virginia Tech
% Silvia Gazzola, University of Bath
% June, 2018

% This file extends the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Get the previous dimension
k = size(T,2)+1;

if k == 1
  v = Atransp_times_vec(A, U(:,k));
  v = v(:);
else
  v = Atransp_times_vec(A, U(:,k));
  v = v(:);
  for j = 1:k-1
    M(j,k)=V(:,j)'*v;
    v = v - M(j,k)*V(:,j);
  end
end
M(k,k) = norm(v);
v = v / M(k,k);

z = L*v;
u = A_times_vec(A, z);
u = u(:);

for j = 1:k
  T(j,k) = U(:,j)'*u;
  u = u - T(j,k)*U(:,j);
end
T(k+1,k) = norm(u);
u = u / T(k+1,k);

U = [U, u];
V = [V, v];

