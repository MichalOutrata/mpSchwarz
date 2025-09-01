function x = trisol_SprsMtrxDnseRhs(T,b,ChopOptions, IsUpperTriang, T_struct)
%trisol_sparse   Solve triangualr system in low precision arithmetic, based on "TRISOL(T, b)", see https://github.com/SrikaraPranesh/Multi_precision_NLA_kernels.
%   TRISOL(T,b,IsUpperTriang) is the solution x to the triangular system
%   Tx = b computed in low precision arithmetic.
%   By default the low precision format is fp16.
%   It requires the "chop" and "chop_sparse" functions to simulate lower precision arithmetic.

if nargin == 2
    fp.format = 'h'; chop_sparse([],fp); chop_dense([],fp);
    if ~norm(T-triu(T),1), IsUpperTriang = true; else, IsUpperTriang = false; end
    if IsUpperTriang
        T_struct = cell(size(T,2),1);
        for col = 1:size(T,2)
            T_struct{col} = find( T(1:col-1,col) );
        end
    else
        T_struct = cell(size(T,2),1);
        for col = 1:size(T,2)
            T_struct{col} = find( T(col+1:end,col) );
        end
    end

elseif nargin == 3
    chop_sparse([],ChopOptions); chop_dense([],ChopOptions);
    if ~norm(T-triu(T),1), IsUpperTriang = true; else, IsUpperTriang = false; end
    if IsUpperTriang
        T_struct = cell(size(T,2),1);
        for col = 1:size(T,2)
            T_struct{col} = find( T(1:col-1,col) );
        end
    else
        T_struct = cell(size(T,2),1);
        for col = 1:size(T,2)
            T_struct{col} = find( T(col+1:end,col) );
        end
    end

elseif nargin == 4
    chop_sparse([],ChopOptions); chop_dense([],ChopOptions);
    if IsUpperTriang
        T_struct = cell(size(T,2),1);
        for col = 1:size(T,2)
            T_struct{col} = find( T(1:col-1,col) );
        end
    else
        T_struct = cell(size(T,2),1);
        for col = 1:size(T,2)
            T_struct{col} = find( T(col+1:end,col) );
        end
    end

else
    chop_sparse([],ChopOptions); chop_dense([],ChopOptions);

end

% Convert matrix and vector to required precisions
b = chop_dense(b); %%% assume that the rhs is NOT sparse -> practical case for DDM subdomain solves
% T = chop_sparse(T); %%% assume that the triangular matrix T is already in precision 'format.fp'


n = length(T); x = zeros(n,1);

if IsUpperTriang        % Upper triangular
    x(n) = chop_dense( b(n)/T(n,n) );
    for i=n-1:-1:1
        temp = chop_sparse( x(i+1) .* T(1:i,i+1));  %%% <- use chop_sparse bcs the most recent column T(1:i,i+1) is sparse (as is T)
        % inds = find(temp); % = T_struct{i+1}                       %%% <- work only on the nnz indcs (ideally we would have rowtree of T and didnt have to call "find")
        b(T_struct{i+1}) = chop_dense( b(T_struct{i+1}) - temp(T_struct{i+1}) ); %%% <- use chop_dense bcs the rhs is assumed dense, ie, the vector is dense
        x(i) = chop_dense( b(i)/T(i,i) );           %%% <- use chop_dense bcs it's just the scalar
        %disp(i + "  " + any(isnan(x)) + "  " + any(isinf(x)))
    end

else                    % Lower triangular
    x(1) = chop_dense( b(1)/T(1,1) );
    for i=2:n
        temp = chop_sparse( x(i-1) .* T(i:n,i-1));  %%% <- use chop_sparse bcs the most recent column T(i:n,i-1) is sparse (as is T)
        % inds = find(temp); % = T_struct{i-1}                       %%% <- work only on the nnz indcs (ideally we would have rowtree of T and didnt have to call "find")
        b(i-1+T_struct{i-1}) = chop_dense( b(i-1+T_struct{i-1}) - temp(T_struct{i-1}) ); %%% <- use chop_dense bcs the rhs is assumed dense, ie, the vector is dense
        x(i) = chop_dense( b(i)/T(i,i) );           %%% <- use chop_dense bcs it's just the scalar
    end


end
