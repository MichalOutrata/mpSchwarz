function [L,U,p,timings_matlab,timings_chop] = lu_sparseMmtrx_chop_timingtester(A,pivoting,format)
%   Sparse triangular factorization **based** on: 
%       textbook [Algorithms for Sparse Linear Systems, Algorithm 6.1; Scott & Tuma, 2023], 
%       paper [A Stable Method for the $LU$ Factorization of M-Matrices; Ahac & Olesky, 1986]
%   For an M-matrix A [L,U,p] = lu_sparseMmtrx_chop(A) produces 
%       a unit lower triangular M-matrix L, 
%       an upper triangular M-matrix U with positive diagonal entries,
%       permutation vector p
%   so that L*U = A(p,p). The permutation vector p is chosen so as to preserve the M-matrix property in a stable way,
%   and that is a dynamical proces, ie, the symbolic factorization cannot be precomputed but rather is computed on the fly. 
%   Moreover, the permutation is generally ignorant of the resulting sparsitry pattern of both L and U.
%
%   It requires the CHOP function to simulate lower precision arithmetic.
%   By default half precision is used

if nargin<2, format = 'h'; end
fp.format = format; chop_sparse([],fp); chop_dense([],fp);

[m,n] = size(A); p = (1:m)';
if strcmp(pivoting,'Mmtrx')
    SumsForCDD = sum(A,1); NeedToPivotForStability = true;
end

timings_matlab = zeros(4,1);
timings_chop = zeros(5,1);

L_initsprs = sparse(size(A));
U_initsprs = sparse(size(A));


for k = 1:n-1

    % find the pivot row/col
    switch(pivoting)
        case 'partial'
          [~,PivotInd_WithinUnrducdMtrx] = max(abs(A(k:n,k))); PivotInd = PivotInd_WithinUnrducdMtrx + k-1;

        case 'Mmtrx'

            if NeedToPivotForStability %%% pivot for stability first and only second for sparsity of the factors
                minCDDpivot = min(SumsForCDD);
                if minCDDpivot < 0 
                    %%% find the new pivot (1st - find positive col-sum columns, 2nd find among them one with lowst #nnzs, 3rd find among them the column with largest pivot)
                    IndsPosCDD = find(SumsForCDD > 0);
                    NnzCols = sum( A(:,IndsPosCDD)~=0 );
                    minNnzCols = min(NnzCols); inds = find(NnzCols == minNnzCols); [~,pivind] = max(diag(A(inds,inds))); PivotInd_WithinUnrducdMtrx = IndsPosCDD(pivind);
                    PivotInd = PivotInd_WithinUnrducdMtrx + k-1;                
                else
                    NeedToPivotForStability = false; 
                    NnzCols = sum( A~=0 ); 
                    minNnzCols = min(NnzCols); inds = find(NnzCols == minNnzCols); [~,PivotInd_WithinUnrducdMtrx] = max(diag(A(inds,inds)));
                    PivotInd = PivotInd_WithinUnrducdMtrx + k-1;
                end

            else %%% pivot for sparsity of the factors only (so called symmetric Markowitz permutation - minimize #nnzs in the pivot column)
                NnzCols = sum( A~=0 ); 
                minNnzCols = min(NnzCols); inds = find(NnzCols == minNnzCols); [~,PivotInd_WithinUnrducdMtrx] = max(diag(A(inds,inds)));
                PivotInd = PivotInd_WithinUnrducdMtrx + k-1;
            end

        case 'none' 
    end

   

    % Skip elimination if column is zero
    switch(pivoting)
                case 'partial'
                    if (A(PivotInd,k) ~= 0)
                        if (PivotInd ~= k)
                            p([k PivotInd]) = p([PivotInd k]);
                            A([k PivotInd],:) = A([PivotInd k],:);
                        end
                    end
                case 'Mmtrx'
                    if (A(PivotInd,k) ~= 0)
                        if (PivotInd ~= k)
                            A([k PivotInd],:) = A([PivotInd k],:); A(:,[k PivotInd]) = A(:,[PivotInd k]);
                            % update the column sums
                            SumsForCDD([1 PivotInd_WithinUnrducdMtrx]) = SumsForCDD([PivotInd_WithinUnrducdMtrx 1]); % update also the column sums "ordering" (swapping rows makes no difference and swapping columns only swaps the two column sums)
                            % SumsForCDD(2:end) = SumsForCDD(2:end) + SumsForCDD(1)*A(1,2:end)/A(1,1); <- this is the formula in [Thm4, Alah-Oleski] but there the off-diag part of mtrx is -a_ij (see (6-7)), so in the normal notation, we need to account for the minus 
                            SumsForCDD(2:end) = SumsForCDD(2:end) + SumsForCDD(1)*(-1)*A(k,k+1:end)/A(k,k); 
                            SumsForCDD = SumsForCDD(2:end);
                        end
                    end
                case 'none' 
    end
    


    % Compute multipliers
    i = k+1:m;
    
    tic
    vec = A(i,k)/A(k,k);
    timings_matlab(1) = timings_matlab(1) + toc;
    
    tic
    vec1 = chop_sparse(vec);
    timings_chop(1) = timings_chop(1) + toc;

    tic
    A(i,k) = vec1;
    timings_matlab(1) = timings_matlab(1) + toc;
    
    % Update the remainder of the matrix
    j = k+1:n;
    
    tic
    mat1 = A(i,k)*A(k,j);
    timings_matlab(2) = timings_matlab(2) + toc;

    tic
    mat2 = chop_sparse(mat1);
    timings_chop(2) = timings_chop(2) + toc;

    tic
    mat3 = A(i,j) - mat2;
    timings_matlab(3) = timings_matlab(3) + toc;  
    
    tic
    A(i,j) = chop_sparse(mat3);
    timings_chop(3) = timings_chop(3) + toc;

  

    tic
    [iaux,jaux,~] = find(mat2);
    auxxx = A(iaux,jaux) - mat2(iaux,jaux);
    timings_matlab(4) = timings_matlab(4) + toc;  

    tic
    auxxx = chop_sparse(auxxx);
    timings_chop(4) = timings_chop(4) + toc;

    tic
    auxxx = chop_dense(full(auxxx));
    timings_chop(5) = timings_chop(5) + toc;

end

% Separate result
L = tril(A,-1) + eye(m,n);
U = triu(A);U = U(1:n,1:end);