function [A,D1,D2,its] = scale_diag_2side_symm(A,tol,prnt)
%SCALE_DIAG_2SIDE_SYMM Symmetry-preserving two-sided diagonal scaling.
%   [B,D1,D2,its] = scale_diag_2side_symm(A,TOL,PRNT) computes
%   B = D1*A_D2 where the diagonal matrices D1 = inv(diag(d1))
%    and D2 = inv(diag(d2)) are chosen
%   so that the largest elements in absolute value in every row and every
%   column of B is 1.  TOL is a convergence tolerance: default TOL = 1e-4.
%   If A has a zero row or column it will stay zero.
%   Set PRNT = 1 (default 0) to print convergence behaviour.
%   ITS is the number of iterations for convergence.

% Reference
% Knight, Ruiz \& U\ccar, A Symmetry Preserving Algorithm for Matrix Scaling,
% SIMAX 35, 931-955, 2014.

%%% adapted code by Outrata, 2025 -> sparse version

if nargin < 2 || isempty(tol), tol = 1e-4; end
if nargin < 3, prnt = 0; end

n = length(A);


d1prod = ones(n,1);
d2prod = ones(1,n);

if prnt
   fprintf('%2.0f:  dA = 0,         normA = %9.2e, condA = %9.2e\n',...
           0, norm(A,1), cond(A,1))
end   

for k = 1:10

    row_max = max(abs(A),[],2);
    col_max = max(abs(A),[],1);
    d1 = sqrt(row_max);
    d2 = sqrt(col_max);
    A_old = A;
    %%% Fixed by Michal
    A = spdiags(1./d1,0,n,n) * A * spdiags(1./d2',0,n,n); %originally missing the transpose, ie, "spdiags(1./d2,0,n,n)" 
    %%% Fixed by Michal
    d1prod = d1prod./d1;
    d2prod = d2prod./d2;
    if prnt
       fprintf('%2.0f:  dA = %9.2e, normA = %9.2e, condA = %9.2e\n',...
               k, norm(A-A_old,1), norm(A,1), cond(A,1))
    end   
    
    row_max = max(abs(A),[],2);
    col_max = max(abs(A),[],1);
    if norm([row_max' col_max] - 1, inf) < tol, its = k; break, end
    
end

D1 = spdiags(d1prod,0,n,n);
%%% Fixed by Michal
D2 = spdiags(d2prod',0,n,n); %originally missing the transpose, ie, "spdiags(1./d2,0,n,n)" 
%%% Fixed by Michal
its = k;