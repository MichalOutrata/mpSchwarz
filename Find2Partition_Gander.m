function [ni,nj,n,m]=Find2Partition_Gander(A)
% FIND2PARRTITION finds a partition of a banded matrix into two blocks
%   [ni,nj,n,m]=Find2Partition(A); finds the ni and nj indices for a
%   partition of the banded matrix A into two submatrices, used for
%   the construction of algebraic optimized Schwarz methods. The
%   partition is then organized as follows:
%
%     nj1 nj2 nj3 nj4  the blocks ni2xnj3 and ni3xnj2 are such that
%      -------------   the banded matrix decouples completely
%  ni1 |  |  |  |  |   
%      -------------   
%  ni2 |  |  |  |  |     
%      ------------- 
%  ni3 |  |  |  |  |
%      -------------
%  ni4 |  |  |  |  |
%      -------------
%   we also return n the size of the matrix, and m the point where
%   it is split

n=size(A,1); m=round(n/2);

i=1;                                % first estimate width of coupling
while isempty(find(A(m+1:end,i), 1))
  i=i+1;      
end
i1=i-1;
i=n;
while isempty(find(A(1:m,i), 1))
  i=i-1;      
end
i2=i+1;
i=1;                                % now estimate height of coupling
while isempty(find(A(i,m+1:end), 1))
  i=i+1;      
end
i3=i-1;
i=n;
while isempty(find(A(i,1:m), 1))
  i=i-1;      
end
i4=i+1;

ni(1)=i3;
ni(2)=m-i3;
ni(3)=i4-m-1;
ni(4)=n-i4+1;

nj(1)=i1;
nj(2)=m-i1;
nj(3)=i2-m-1;
nj(4)=n-i2+1;

end
