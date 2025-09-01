function [x,timings_chop,timings_chopsprs] = trisol_sparse_timingtester(T,b,IsUpperTriang, format)
%trisol_sparse   Solve triangualr system in low precision arithmetic, based on "TRISOL(T, b)", see https://github.com/SrikaraPranesh/Multi_precision_NLA_kernels.
%   TRISOL(T,b,IsUpperTriang) is the solution x to the triangular system
%   Tx = b computed in low precision arithmetic.
%   By default the low precision format is fp16.
%   It requires the "chop" and "chop_sparse" functions to simulate lower precision arithmetic.

if nargin < 4
    fp.format = 'h'; chop_sparse([],fp); chop_dense([],fp);
else
    fp.format = format; chop_sparse([],fp); chop_dense([],fp);
end

timings_chop = zeros(3,1);
timings_chopsprs = zeros(4,1);


% Convert matrix and vector to required precisions
tic
b = chop_dense(b); %%% assume that the rhs is NOT sparse -> practical case for DDM subdomain solves
timings_chop(1) = toc;
tic
T = chop_sparse(T); %%% assume that the triangular matrix T is already in precision 'format.fp'
timings_chopsprs(1) = toc;

n = length(T); x = zeros(n,1);

if IsUpperTriang        % Upper triangular
   
    tic
    x(n) = chop_dense( b(n)/T(n,n) );
    timings_chop(2) = timings_chop(2) + toc;

    for i=n-1:-1:1

        tic
        temp = chop_sparse( x(i+1) .* T(1:i,i+1));  %%% <- use chop_sparse bcs the most recent column T(1:i,i+1) is sparse (as is T)
        timings_chopsprs(2) = timings_chopsprs(2) + toc;
        
        tic
        b(1:i) = chop_dense( b(1:i) - temp );       %%% <- use chop_dense bcs the rhs is assumed dense, ie, the vector is dense
        timings_chop(3) = timings_chop(3) + toc;

        tic
        x(i) = chop_dense( b(i)/T(i,i) );           %%% <- use chop_dense bcs it's just the scalar
        timings_chop(2) = timings_chop(2) + toc;
    end



else                    % Lower triangular
    tic
    x(1) = chop_dense( b(1)/T(1,1) );
    timings_chop(2) = timings_chop(2) + toc;

    for i=2:n
        
        tic
        temp = chop_sparse( x(i-1) .* T(i:n,i-1));  %%% <- use chop_sparse bcs the most recent column T(i:n,i-1) is sparse (as is T)
        timings_chopsprs(2) = timings_chopsprs(2) + toc;
        
        tic
        b(i:n) = chop_dense( b(i:n) - temp );       %%% <- use chop_dense bcs the rhs is assumed dense, ie, the vector is dense
        timings_chop(3) = timings_chop(3) + toc;

        tic
        inds = find(temp);
        b(inds) = chop_dense( b(inds) - temp(inds) );       %%% <- use chop_dense bcs the rhs is assumed dense, ie, the vector is dense
        timings_chopsprs(3) = timings_chopsprs(3) + toc;

        tic
        b(inds) = chop_dense( b(inds) - temp(inds) );       %%% <- use chop_dense bcs the rhs is assumed dense, ie, the vector is dense
        timings_chopsprs(4) = timings_chopsprs(4) + toc;


        tic
        x(i) = chop_dense( b(i)/T(i,i) );           %%% <- use chop_dense bcs it's just the scalar
        timings_chop(2) = timings_chop(2) + toc;
    end


end
