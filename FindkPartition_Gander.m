function [ni,nj,n,m]=FindkPartition_Gander(A,k)
% FINDKPARRTITION finds a partition of a banded matrix into 2^k blocks
%   [ni,nj,n,m]=FindkPartition(A); finds the ni and nj indices for a
%   partition of the banded matrix A into 2^k submatrices (recursively),
%   using Find2Partition().

[ni,nj,n,m] = Find2Partition_Gander(A);
if(k>1)
    [ni1,nj1,n1,m1] = FindkPartition_Gander(A(1:m,1:m),k-1);
    [ni2,nj2,n2,m2] = FindkPartition_Gander(A(m+1:end,m+1:end),k-1);
    ni = [ni1(1:end-1) ni1(end)-ni(2) ni(2) ni(3) ni2(1)-ni(3) ni2(2:end)];
    nj = [nj1(1:end-1) nj1(end)-nj(2) nj(2) nj(3) nj2(1)-nj(3) nj2(2:end)];
end
