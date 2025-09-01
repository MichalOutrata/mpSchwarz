% B = [-1,2;3,4];
% foo = isDiagDom_test( B )
% B = [-1,2;0.5,4];
% foo = isDiagDom_test( B )
% B = [-1,2;-1.01,4];
% foo = isDiagDom_test( B,1e-1 )



function [ AisDiagDom ] = isDiagDom( A, varargin )

AisDiagDom = false; [~,n] = size(A);
if nargin == 1, tol = 1e-15; else, tol = varargin{1}; end
for col = 1:n
    diag_abs_val = abs(A(col,col)); offdiag_absval_sum = sum(abs(A(1:col-1,col)))+sum(abs(A(col+1:n,col)));
    if diag_abs_val < offdiag_absval_sum  &&  abs(diag_abs_val - offdiag_absval_sum) > tol
        col = col-1; break; end
end
if col == n, AisDiagDom = true; end

end

