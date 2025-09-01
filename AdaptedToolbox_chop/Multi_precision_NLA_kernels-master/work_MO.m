%WORK1 tester for various multiprecision kernels
clear all; close all; clc

n = 5;
A = diag(ones(n,1),0) - diag(0.5*ones(n-1,1),-1) - ( [ones(n-1,1);0] * [zeros(n-1,1)', 1] ); Ainv = inv(A);
[l,u] = lu(A);

perm = [5 1 2 3 4]; A_perm = A(perm,perm);

[L,U,p] = lu_sparseMmtrx_chop(A_perm,'Mmtrx');


clear all; close all; clc

n = 4000; nmb_inds = 200;

den = 10/n; a = sprand(n,n,den);
inds = randi([1 n],1,nmb_inds);
inds_row = randi([1 n],1,nmb_inds); 
inds_col = randi([1 n],1,nmb_inds); 

tic
aux = a(inds,inds);
time_randinds = toc

tic
aux = a(min(inds):min(inds)+nmb_inds,min(inds)+10:min(inds)+20+nmb_inds);
time_arrayinds = toc

tic
aux = a(inds_row,inds_col) - 2*a(inds_row,inds_col);
time_substrinds = toc


clear all; close all; clc

n = 320; a = delsq(numgrid('S',n));
[l,u,p,q] = lu(a);

tic
[count,h,parent,post,L] = symbfact(l,'lo','lower');
time_symbfact = toc