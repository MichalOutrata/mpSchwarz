% B = [-1,2;3,4];
% foo = isMmtrx_test( B )


function [ AisMmtrx ] = isMmtrx( A )

% condition 1 - positive diagonal
AisMmtrx = false;
d = diag(A); if min(d) < 0, AisMmtrx = false; return; end

% condition 2 - non-positive off-diagonal
[~,~,A_offdiag_nnzs] = find( spdiags(zeros(size(A,1),1),0,A) );
if max(A_offdiag_nnzs) > 0, AisMmtrx = false; return; end

%%% AisZmtrx = true

% condition 3 - every leading principal minor is positive (Bermmann&Plemmons, Chapt.6, Thm 2.3.A1)
% for ind = 1:size(A,1)
%     if det(A(1:ind,1:ind)) < 0, AisMmtrx = false; return; end;
% end

% condition 3 - the spectral radius of B := I - diag(A)^{-1}*A is smaller than one, i.e., rho( B ) < 1 
% tic
B_aux = spdiags(zeros(size(A,1),1),0, diag(diag(A))\A); eigmax = eigs(B_aux,1,'largestabs','MaxIterations',1000);
if abs(eigmax) < 1, AisMmtrx = true; end % disp('AisMmtrx \w eigval test = true'); end
% toc

% tic
% for ind = 1:size(A,1)
%     if det(A(1:ind,1:ind)) < 0, AisMmtrx = false; return; end;
% end
% toc
% disp('AisMmtrx \w principal leading minors = true');


end

