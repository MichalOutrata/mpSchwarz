function[MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethod,RoundingStuff, debug)

SM_types = SchwarzMethod{1}; 
SM_nmbsubdom_input = SchwarzMethod{2}; 
CalcErrMtrcs = false; %%% change if we want to check some theoretical conds also while running mpSM as preconditioners


%%% how many subdomains
N = size(A,1);
SM_nmbsubdom_PwrOfTwo = floor(log2(SM_nmbsubdom_input)); SM_nmbsubdom = 2^SM_nmbsubdom_PwrOfTwo;
if SM_nmbsubdom ~= SM_nmbsubdom_input, disp(append(append(append('The nmb of subdoms needs to be 2^k for some k. Hence we changed the given',num2str(SM_nmbsubdom_input)),' to '),num2str(SM_nmbsubdom))); end


%%% do the partitioning
[ni,nj,~,~] = FindkPartition_Gander(A,SM_nmbsubdom_PwrOfTwo); rows_cumsum = cumsum(ni); cols_cumsum = cumsum(nj);
frst_rows = NaN(SM_nmbsubdom,1); last_rows = NaN(SM_nmbsubdom,1);  frst_cols = NaN(SM_nmbsubdom,1); last_cols = NaN(SM_nmbsubdom,1);
DoWeRun_RAS = false;
for ind = 1:length(SM_types)
    if strcmp(SM_types{ind},'RAS'), mid_rows = NaN(SM_nmbsubdom,1); mid_cols = NaN(SM_nmbsubdom,1); DoWeRun_RAS = true; end
end

for ind_subdom = 1:SM_nmbsubdom
    if ind_subdom == 1
        frst_rows(1) = 1; frst_cols(1) = 1;
    else
        frst_rows(ind_subdom) = rows_cumsum( (ind_subdom-1-1)*3 +1 ) + 1; frst_cols(ind_subdom) = cols_cumsum( (ind_subdom-1-1)*3 +1 ) + 1;
    end

    if DoWeRun_RAS % if we do RAS -> get the midpoints of the subdomain overlaps
        if ind_subdom ~= SM_nmbsubdom
            mid_rows(ind_subdom) = rows_cumsum( (ind_subdom-1)*3 +2 ); mid_cols(ind_subdom) = cols_cumsum( (ind_subdom-1)*3 +2 );
        end
    end
 
    if ind_subdom ~= SM_nmbsubdom
        last_rows(ind_subdom) = rows_cumsum( (ind_subdom-1)*3 +3 ); last_cols(ind_subdom) = cols_cumsum( (ind_subdom-1)*3 +3 );
    else
        last_rows(ind_subdom) = N; last_cols(ind_subdom) = N;
    end
end

if DoWeRun_RAS % if we do RAS -> get the restriction indices
    frst_rows_RAS = NaN(SM_nmbsubdom,1); last_rows_RAS = NaN(SM_nmbsubdom,1);
    for ind_subdom = 1:SM_nmbsubdom
        if ind_subdom == 1
            frst_rows_RAS(ind_subdom) = 1; last_rows_RAS(ind_subdom) = mid_rows(ind_subdom);
        elseif ind_subdom == SM_nmbsubdom
            frst_rows_RAS(ind_subdom) = mid_rows(ind_subdom-1)+1; last_rows_RAS(ind_subdom) = N;
        else
            frst_rows_RAS(ind_subdom) = mid_rows(ind_subdom-1)+1; last_rows_RAS(ind_subdom) = mid_rows(ind_subdom);
        end
    end
end


%%% do the LowPrecision blocks
MySubdomFactors_LP = cell(SM_nmbsubdom,1); % need to keep the L,U,P,Q as well as other things (eg, A_diag for Stieltjess-type rounding)

for ind_subdom = 1:SM_nmbsubdom
    MySubdomFactors_LP{ind_subdom} = SubdomProbRounding_ChopOrAdvnpx( A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom)) , RoundingStuff, CalcErrMtrcs );
end


%%% pacakage the output
if DoWeRun_RAS, MySubdomIndices = {frst_rows,last_rows,frst_rows_RAS,last_rows_RAS}; else, MySubdomIndices = {frst_rows,last_rows}; end
MyPrecondPackage = {SchwarzMethod,RoundingStuff,MySubdomFactors_LP,MySubdomIndices};

end