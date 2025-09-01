function x = FrwrdBckwrdSubst_MPSM(StoredFactors,b,format)
%   Carry out forward and backward substitution for the subdomain systems in low precision arithmetic.
%   Based on sparse adaptation of "TRISOL(T, b)", see https://github.com/SrikaraPranesh/Multi_precision_NLA_kernels
%   By default the low precision format is fp16.
%   It requires the "chop_dense" and "chop_sparse" functions to simulate lower precision arithmetic.

if nargin == 2
    fp.format = 'h'; chop_sparse([],fp); chop_dense([],fp);
else
    fp.format = format; chop_sparse([],fp); chop_dense([],fp);

end

% Convert rhs to required precisions (assume that the rhs is NOT sparse -> practical case for DDM subdomain solves)
mu_FitToRange_Rhs = StoredFactors{1} / max(abs(b));
b = chop_dense(mu_FitToRange_Rhs * b);  b_bfrperm = b_bfrperm; n = length(b);

if length(StoredFactors) == 10 %%% We deal with the M-matrix rounding case
    

    %%% Solve with "L"
    y = zeros(n,1);
    y(1) = chop_dense( b(1)/L(1,1) );
    for i=2:n
        temp = chop_sparse( y(i-1) .* L(i:n,i-1));  %%% <- use chop_sparse bcs the most recent column L(i:n,i-1) is sparse (as is L)
        b(L_struct{i-1}) = chop_dense( b(L_struct{i-1}) - temp(L_struct{i-1}) ); %%% <- use chop_dense bcs the rhs is assumed dense, ie, the vector is dense
        y(i) = chop_dense( b(i)/L(i,i) );           %%% <- use chop_dense bcs it's just the scalar
    end

    %%% Solve with "U"
    x = zeros(n,1);
    x(n) = chop_dense( y(n)/U(n,n) );
    for i=n-1:-1:1
        temp = chop_sparse( x(i+1) .* U(1:i,i+1));  %%% <- use chop_sparse bcs the most recent column U(1:i,i+1) is sparse (as is U)
        y(U_struct{i+1}) = chop_dense( y(U_struct{i+1}) - temp(T_struct{i+1}) ); %%% <- use chop_dense bcs the rhs is assumed dense, ie, the vector is dense
        x(i) = chop_dense( y(i)/U(i,i) );           %%% <- use chop_dense bcs it's just the scalar
    end
    

    
    









elseif length(StoredFactors) == 11 %%% We deal with the SandwichScaling case








end




