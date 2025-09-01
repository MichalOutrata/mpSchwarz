% B = [-1,2;3,4];
% foo = isSymm_test( B )
% B = [-1,2;2,4];
% foo = isSymm_test( B )
% B = [-1,1e-4;2*1e-4,4];
% foo = isSymm_test( B )
% B = [-1,1+1e-4;1-2*1e-4,4];
% foo = isSymm_test( B, 1e-3 )
% B = [-1,1+1e-4;1-2*1e-4,4];
% foo = isSymm_test( B, 1e-4 )


function [ AisSymm ] = isSymm( A, varargin )

AisSymm = false;
if ~isreal(A), disp('Input is not a real matrix. We will check A=A^T, not A=A^*. Continue?'); waitforbuttonpress; end
if nargin == 1, tol = 1e-14; else, tol = varargin{1}; end
[~,~,At_vec] = find(A'); [~,~,A_vec] = find(A); E_vec = A_vec - At_vec;
if all( abs(E_vec)./abs(A_vec) < tol ), AisSymm = true; end

end

