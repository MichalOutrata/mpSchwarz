% clear; clc; close('all');
% 
% %%% choose set-up
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Schwarz set-up
% SM_type = "RAS"; % "RAS"; % "MS"; % which methods we want to run out of "dAS", "RAS", "MS" 
% SM_nmbsubdoms_PwrOfTwo = 1; % number of subdomains for the Schwarz method
% SM_nmbiter = 1; % number of RAS iterations
% SM_relresacc = 1e-12; % number of RAS iterations
% dampingTheta = 1/3; % generally one should use 1/(2^SM_nmbsubdoms_PwrOfTwo+1) but since we have sausage-like domain decomposition, we can do with alternating coloring no matter how many subdomains
% SM_RandomRHS = true;
% 
% %%% RoundingStuff set-up
% RM_Software = "advanpix";   % "advanpix" or "chop"
% RM_Type = "Mmtrx";          % "Mmtrx" or "StieltjessDiagDom"
% if strcmp(RM_Software,"advanpix")
%     RM_PrecisionChoice = 5; % number of digits advanpix uses for the subdomain solves
% elseif strcmp(RM_Software,"chop")
%     ChopOptions.explim = 'fp16'; % precision format of chop for the subdomain solves  
%     ChopOptions.explim = 1; % do we consider the "range limitation". If "0" then bound on the maximal exponent is ignored inside chop   
%     RM_PrecisionChoice = ChopOptions;
% end
% RM_ScalingToRange = true;           % whether or not we rescaled the matrix to fit to the considered low-precision range of representation 
% RM_RescaleType = "diag2side";       % "diag2side" or "diag2side_symm"
% RM_ScalingToRange_Fraction = 16;    % we rescale the matrix and RHS to fit to "x_max / RM_ScalingToRange_Fraction"
% RM_ScalingOptions = {RM_ScalingToRange,RM_RescaleType,RM_ScalingToRange_Fraction};
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nmb_int_gridcols = 50; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
% h=1/(nmb_int_gridcols+1); % mesh size
% 
% %%% do non-symmetric AdvecDiff based on SzyldGander Fig 2.1
% G=numgrid('S',nmb_int_gridcols+2);
% eta = @(x,y) x.^2.*cos(x+y).^2'; a = @(x,y) 20*(x+y).^2.*exp(x-y); b1 =  @(x,y) (y-0.5); b2 =  @(x,y) -(x-0.5);
% A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,a,b1,b2});
% u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; x_mesh = 0:h:1; 
% v = ones(size(A,1),1);
% debug_bool = false;
% 
% %%% run test
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % We test whether one iteration of a given SM for "u_init = 0" & "b" IS THE SAME as application of the given SM preconditioner   
% SM_type = "dAS";
% if strcmp(SM_type,'dAS'), SchwarzMethod = {SM_type, 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc, dampingTheta}; else, SchwarzMethod = {SM_type, 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc}; end
% RM_RoundingOptions = {RM_PrecisionChoice,RM_ScalingOptions}; RoundingStuff = {RM_Software, RM_Type, RM_RoundingOptions};
% 
% [MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethod, RoundingStuff, debug_bool);
% [t,precvec] = mpSchwarzMethodsPreconds_test(v, {A, MyPrecondPackage});
% [MyOutput] = mpSchwarzMethods(A,v,zeros(size(v)), SchwarzMethod, RoundingStuff, {false,false,false}, {}, debug_bool);
% str = "Pro " + SM_type + " ~ " + norm(precvec - MyOutput{1}')/norm(MyOutput{1}); disp(str)
% 
% SM_type = "RAS";
% if strcmp(SM_type,'dAS'), SchwarzMethod = {SM_type, 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc, dampingTheta}; else, SchwarzMethod = {SM_type, 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc}; end
% RM_RoundingOptions = {RM_PrecisionChoice,RM_ScalingOptions}; RoundingStuff = {RM_Software, RM_Type, RM_RoundingOptions};
% 
% [MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethod, RoundingStuff, debug_bool);
% [t,precvec] = mpSchwarzMethodsPreconds_test(v, {A, MyPrecondPackage});
% [MyOutput] = mpSchwarzMethods(A,v,zeros(size(v)), SchwarzMethod, RoundingStuff, {false,false,false}, {}, debug_bool);
% str = "Pro " + SM_type + " ~ " + norm(precvec - MyOutput{1}')/norm(MyOutput{1}); disp(str)
% 
% SM_type = "MS";
% if strcmp(SM_type,'dAS'), SchwarzMethod = {SM_type, 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc, dampingTheta}; else, SchwarzMethod = {SM_type, 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc}; end
% RM_RoundingOptions = {RM_PrecisionChoice,RM_ScalingOptions}; RoundingStuff = {RM_Software, RM_Type, RM_RoundingOptions};
% 
% [MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethod, RoundingStuff, debug_bool);
% [t,precvec] = mpSchwarzMethodsPreconds_test(v, {A, MyPrecondPackage});
% [MyOutput] = mpSchwarzMethods(A,v,zeros(size(v)), SchwarzMethod, RoundingStuff, {false,false,false}, {}, debug_bool);
% str = "Pro " + SM_type + " ~ " + norm(precvec - MyOutput{1}')/norm(MyOutput{1}); disp(str)
% 
% 






function[t_out,Prec_v] = mpSchwarzMethodsPreconds(v, AdditionalArgs)

A = AdditionalArgs{1}; mpSchwarzPrecond = AdditionalArgs{2};
SchwarzMethod = mpSchwarzPrecond{1}; RoundingStuff = mpSchwarzPrecond{2}; MySubdomFactors_LP = mpSchwarzPrecond{3}; MySubdomIndices = mpSchwarzPrecond{4};
SM_type = SchwarzMethod{1}; SM_nmbsubdom_input = SchwarzMethod{2}; if strcmp(SM_type,'dAS'), damping_theta = SchwarzMethod{end}; end %; SM_nmbiter = SchwarzMethod{3}; SM_relresacc = SchwarzMethod{4}
% we have " RoundingStuff = {RoundingSoftware, RoundingType, RoundingOptions} "
RoundingSoftware = RoundingStuff{1};    %%% "advanpix" or "chop"
RoundingType = RoundingStuff{2};        %%% type of rounding routine we call - "Mmtrx" or "StieltjessDiagDom"
RoundingOptions = RoundingStuff{3};     %%% extra params for rounding software (notably ChopOptions/NmbDigits for chop/advanpix) 
ScalingOptions = RoundingOptions{2};    %%% additional parameters for the rounding software and/or the functions
RescaleToFitRange = ScalingOptions{1};  %%% whether or not we rescaled the matrix to fit to the considered low-precision range of representation  


SM_nmbsubdom_PwrOfTwo = floor(log2(SM_nmbsubdom_input)); SM_nmbsubdom = 2^SM_nmbsubdom_PwrOfTwo; 
if SM_nmbsubdom ~= SM_nmbsubdom_input, disp(append(append(append('The nmb of subdoms needs to be 2^k for some k. Hence we changed the given',num2str(SM_nmbsubdom_input)),' to '),num2str(SM_nmbsubdom))); end
N = length(v); Prec_v = zeros(N,1);
frst_rows = MySubdomIndices{1}; last_rows = MySubdomIndices{2}; 
if strcmp(SM_type,'RAS'), frst_rows_RAS = MySubdomIndices{3}; last_rows_RAS = MySubdomIndices{4}; end
res_prev = v;

TimerSM = tic;
for ind_subdom = 1:SM_nmbsubdom
    
    %%% calculate the subdomain solution  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    curr_rhs_bfrPerm_bfrScale = res_prev(frst_rows(ind_subdom):last_rows(ind_subdom));        
    
    if strcmp(RoundingSoftware,'chop') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % "RescaleToFitRange" - if we rescaled the mtrx bfr factorization 
        % "RoundingType" - type of rounding routine we call - "Mmtrx" or "StieltjessDiagDom" 
        % my current output of the low-prec factorization is stored in "MySubdomFactors_LP{ind_subdom}{:}"
        % ChopOutput = {mu_FitToRange, D_rowscale,D_colscale, L,L_struct, U,U_struct, p,q_inv, ErrMtrx,A_Rounded,A_RescaledToFitRange}; 
        ChopOptions = RoundingOptions{1}; chop_dense([],ChopOptions);
    
        if RescaleToFitRange
            curr_mu_FitToRange = MySubdomFactors_LP{ind_subdom}{1}; curr_D_rowscale = MySubdomFactors_LP{ind_subdom}{2}; curr_D_colscale = MySubdomFactors_LP{ind_subdom}{3};
            % scale the rhs vector
            curr_rhs_bfrPerm_SemiScaled = curr_D_rowscale * curr_rhs_bfrPerm_bfrScale;
            if strcmp('fp64',ChopOptions.format) || strcmp('fp32',ChopOptions.format) || strcmp('fp16',ChopOptions.format)
                ExtraFactorForRHS = 100;
            elseif strcmp('bfloat16',ChopOptions.format)
                ExtraFactorForRHS = 10;
            elseif strcmp('q43',ChopOptions.format) || strcmp('q52',ChopOptions.format)
                ExtraFactorForRHS = 5;
            end
            curr_rhsScalingFactor = max(curr_rhs_bfrPerm_SemiScaled)*ExtraFactorForRHS;
            curr_rhs_bfrPerm_Scaled = (curr_rhs_bfrPerm_SemiScaled/curr_rhsScalingFactor) * curr_mu_FitToRange;
            % permute the rhs vector
            curr_Pvec = MySubdomFactors_LP{ind_subdom}{8}; curr_rhs_Permed_Scaled = curr_rhs_bfrPerm_Scaled(curr_Pvec);
            % solve with the precomputed factors
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{4}; curr_Lfact_NnzStruct = MySubdomFactors_LP{ind_subdom}{5}; 
            curr_Ufact = MySubdomFactors_LP{ind_subdom}{6}; curr_Ufact_NnzStruct = MySubdomFactors_LP{ind_subdom}{7}; 
    
            curr_rhs_chopped = chop_dense(curr_rhs_Permed_Scaled); 
            uaux = trisol_SprsMtrxDnseRhs(curr_Lfact,curr_rhs_chopped, ChopOptions,false,curr_Lfact_NnzStruct);
            curr_u_bfrPerm_bfrScale = trisol_SprsMtrxDnseRhs(curr_Ufact,uaux, ChopOptions,true,curr_Ufact_NnzStruct);
            % permute the solution vector back
            curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{9}; curr_u_Permed_bfrScale = double(curr_u_bfrPerm_bfrScale(curr_Qvec_inv));
            curr_u_Permed_Scaled = curr_D_colscale * (curr_u_Permed_bfrScale*curr_rhsScalingFactor);
    
        else
            % permute the rhs vector
            curr_Pvec = MySubdomFactors_LP{ind_subdom}{8}; curr_rhs_Permed = curr_rhs_bfrPerm_bfrScale(curr_Pvec);
            % solve with the precomputed factors
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{4}; curr_Lfact_NnzStruct = MySubdomFactors_LP{ind_subdom}{5}; 
            curr_Ufact = MySubdomFactors_LP{ind_subdom}{6}; curr_Ufact_NnzStruct = MySubdomFactors_LP{ind_subdom}{7}; 
    
            curr_rhs_chopped = chop_dense(curr_rhs_Permed); 
            %norm(curr_rhs_chopped-curr_rhs_Permed)/norm(curr_rhs_Permed)
            uaux = trisol_SprsMtrxDnseRhs(curr_Lfact,curr_rhs_chopped, ChopOptions,false,curr_Lfact_NnzStruct);
            curr_u_bfrPerm_bfrScale = trisol_SprsMtrxDnseRhs(curr_Ufact,uaux, ChopOptions,true,curr_Ufact_NnzStruct);
    
            % permute the solution vector back
            curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{9}; curr_u_Permed = double(curr_u_bfrPerm_bfrScale(curr_Qvec_inv));
            curr_u_Permed_Scaled = curr_u_Permed;
        end
    
    
    elseif strcmp(RoundingSoftware,'advanpix')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % "RescaleToFitRange" - if we rescaled the mtrx bfr factorization 
        % "RoundingType" - type of rounding routine we call - "Mmtrx" or "StieltjessDiagDom" 
        % "Advanpix_NmbDigits = RoundingOptions{1}" - number of digits used for our low-precision 
        % my current output of the low-prec factorization is stored in "MySubdomFactors_LP{ind_subdom}{:}"
        % AdvanpixOutput = {mu_FitToRange, D_rowscale,D_colscale, L, U, p,q_inv, ErrMtrx_mp,A_rounded,A_RescaledToFitRange}; 
    
        if RescaleToFitRange
            curr_mu_FitToRange = MySubdomFactors_LP{ind_subdom}{1}; curr_D_rowscale = MySubdomFactors_LP{ind_subdom}{2}; curr_D_colscale = MySubdomFactors_LP{ind_subdom}{3};
            % scale the rhs vector
            curr_rhs_bfrPerm_SemiScaled = curr_D_rowscale * curr_rhs_bfrPerm_bfrScale;
            curr_rhsScalingFactor = max(curr_rhs_bfrPerm_SemiScaled)*2;
            curr_rhs_bfrPerm_Scaled = (curr_rhs_bfrPerm_SemiScaled/curr_rhsScalingFactor) * curr_mu_FitToRange;
            % permute the rhs vector
            curr_Pvec = MySubdomFactors_LP{ind_subdom}{6}; curr_rhs_Permed_Scaled = curr_rhs_bfrPerm_Scaled(curr_Pvec);
            % solve with the precomputed factors
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{4}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{5}; Advanpix_NmbDigits = RoundingOptions{1};
            warning('off','all') %%% is it really necessary?
            mp.Digits(Advanpix_NmbDigits); curr_rhs_mp = mp(curr_rhs_Permed_Scaled); curr_u_bfrPerm_bfrScale = curr_Ufact \ ( curr_Lfact \ curr_rhs_mp); 
            warning('on','all')
            % permute the solution vector back
            curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{7}; curr_u_Permed_bfrScale = double(curr_u_bfrPerm_bfrScale(curr_Qvec_inv));
            % scale the solution vector back
            curr_u_Permed_Scaled = curr_D_colscale * (curr_u_Permed_bfrScale*curr_rhsScalingFactor);
    
        else
            % permute the rhs vector
            curr_Pvec = MySubdomFactors_LP{ind_subdom}{6}; curr_rhs_Permed = curr_rhs_bfrPerm_bfrScale(curr_Pvec);
            % solve with the precomputed factors
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{4}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{5}; Advanpix_NmbDigits = RoundingOptions{1};
            warning('off','all') %%% is it really necessary?
            mp.Digits(Advanpix_NmbDigits); curr_rhs_mp = mp(curr_rhs_Permed); curr_u_bfrPerm_bfrScale = curr_Ufact \ ( curr_Lfact \ curr_rhs_mp); 
            warning('on','all')
            % permute the solution vector back
            curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{7}; curr_u_Permed = double(curr_u_bfrPerm_bfrScale(curr_Qvec_inv));
            curr_u_Permed_Scaled = curr_u_Permed;
        end

    end





    %%% the subdomain updating strategies, depending on Schwarz method type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ind_subdom == 1 %%% update my global solution for subdomain 1
        if strcmp(SM_type,'RAS') 
            Prec_v = Prec_v + [curr_u_Permed_Scaled(1:last_rows_RAS(ind_subdom)); zeros(N-last_rows_RAS(ind_subdom),1)];
        elseif strcmp(SM_type,'dAS')
            Prec_v = Prec_v + damping_theta * [curr_u_Permed_Scaled; zeros(N-last_rows(ind_subdom),1)];
        elseif strcmp(SM_type,'MS')
            Prec_v = Prec_v + [curr_u_Permed_Scaled; zeros(N-last_rows(ind_subdom),1)];
        end

    elseif ind_subdom == SM_nmbsubdom %%% update my global solution for subdomain N
        if strcmp(SM_type,'RAS') 
            frst_row_inside_subdom = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1;
            Prec_v = Prec_v + [zeros(last_rows_RAS(ind_subdom-1),1); curr_u_Permed_Scaled(frst_row_inside_subdom:end)]; 
        elseif strcmp(SM_type,'dAS')
            Prec_v = Prec_v + damping_theta * [zeros(frst_rows(ind_subdom)-1,1); curr_u_Permed_Scaled]; 
        elseif strcmp(SM_type,'MS')
            Prec_v = Prec_v + [zeros(frst_rows(ind_subdom)-1,1); curr_u_Permed_Scaled]; 
        end

    else %%% update my global solution for middle subdomains
        if strcmp(SM_type,'RAS') 
            frst_row_inside_subdom = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1; % nmb_of_rows_in_left_overlap = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) -> doesn't include the +1 to start after the overlap
            last_row_inside_subdom = last_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1; % nmb_of_rows_to_right_overlap = last_rows_RAS(ind_subdom) - frst_rows(ind_subdom) -> doesn't include the +1 to start after the overlap
            Prec_v = Prec_v + [zeros(last_rows_RAS(ind_subdom-1),1); curr_u_Permed_Scaled(frst_row_inside_subdom:last_row_inside_subdom); zeros(N-last_rows_RAS(ind_subdom),1)];
        elseif strcmp(SM_type,'dAS')
            Prec_v = Prec_v + damping_theta * [zeros(frst_rows(ind_subdom)-1,1); curr_u_Permed_Scaled; zeros(N-last_rows(ind_subdom),1)];
        elseif strcmp(SM_type,'MS')
            Prec_v = Prec_v + [zeros(frst_rows(ind_subdom)-1,1); curr_u_Permed_Scaled; zeros(N-last_rows(ind_subdom),1)];
        end
    end


    if strcmp(SM_type,'MS') && ind_subdom ~= SM_nmbsubdom % if multiplicative Schwarz -> update the residual after each subdomain solve
        res_prev = v-A*Prec_v; 
    end

end
    
t = toc(TimerSM); t_out = [t;0];
end