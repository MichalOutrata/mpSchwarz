function[MyOutput] = mpSchwarzMethods(A,rhs,u_init, SchwarzMethod, RoundingStuff, ErrorCheckStuff, plot,debug)

    function v_out = MatVecFun(v_in,tflag, Msolve,Mmult)
    if strcmp(tflag,'notransp'), v_out_aux = Mmult * v_in; v_out = Msolve \ v_out_aux;
    else, v_out_aux = Msolve' \ v_in; v_out = Mmult' * v_out_aux; end
    end

SM_type = SchwarzMethod{1}; SM_nmbsubdom_input = SchwarzMethod{2}; SM_nmbiter = SchwarzMethod{3}; SM_relresacc = SchwarzMethod{4}; if strcmp(SM_type,'dAS'), damping_theta = SchwarzMethod{end}; end

% we have " RoundingStuff = {RoundingSoftware, RoundingType, RoundingOptions} "
RoundingSoftware = RoundingStuff{1};    %%% "advanpix" or "chop"
RoundingType = RoundingStuff{2};        %%% type of rounding routine we call - "Mmtrx" or "StieltjessDiagDom"
RoundingOptions = RoundingStuff{3};     %%% extra params for rounding software (notably ChopOptions/NmbDigits for chop/advanpix) 
ScalingOptions = RoundingOptions{2};    %%% additional parameters for the rounding software and/or the functions
RescaleToFitRange = ScalingOptions{1};  %%% whether or not we rescaled the matrix to fit to the considered low-precision range of representation  
CheckErrMtrcs_norms = ErrorCheckStuff{1};  %%% whether or not we check the condition ||mathcalAi_inv * F_i || < 1  
CheckErrMtrcs_entrs = ErrorCheckStuff{2};  %%% whether or not we check the condition mathcalAi_inv >= mathcalAi_inv * F_i * mathcalAi_inv   
CheckErrMtrcs_PregSplit = ErrorCheckStuff{3};  %%% whether or not we check the condition mathcalAi_inv >= mathcalAi_inv * F_i * mathcalAi_inv   

if isempty(plot)
    plot_convergence = 0;
else    
    plot_convergence = 1; x_mesh = plot{1}; u_plot = plot{2}; waittime = plot{3}; angle1 = plot{4}; angle2 = plot{5};
end


N = length(rhs); u_seq = NaN(SM_nmbiter,N);
SM_nmbsubdom_PwrOfTwo = floor(log2(SM_nmbsubdom_input)); SM_nmbsubdom = 2^SM_nmbsubdom_PwrOfTwo;
if SM_nmbsubdom ~= SM_nmbsubdom_input, disp(append(append(append('The nmb of subdoms needs to be 2^k for some k. Hence we changed the given',num2str(SM_nmbsubdom_input)),' to '),num2str(SM_nmbsubdom))); end
if CheckErrMtrcs_norms || CheckErrMtrcs_entrs || CheckErrMtrcs_PregSplit, CalcErrMtrcs = true; else, CalcErrMtrcs = false; end
if CheckErrMtrcs_norms, ErrMtrcs_NormCond = false(SM_nmbsubdom,1); ErrMtrcs_NormVals = zeros(SM_nmbsubdom,1); else, ErrMtrcs_NormCond = nan; ErrMtrcs_NormVals = nan; end
if CheckErrMtrcs_entrs, ErrMtrcs_EntrCond = false(SM_nmbsubdom,1); ErrMtrcs_EntrsViolation = zeros(SM_nmbsubdom,1); else, ErrMtrcs_EntrCond = nan; ErrMtrcs_EntrsViolation = nan; end
if CheckErrMtrcs_PregSplit, ErrMtrcs_PregSplitCond = false(SM_nmbsubdom,1); ErrMtrcs_PregSplitEigRatio = zeros(SM_nmbsubdom,1); else, ErrMtrcs_PregSplitCond = nan; ErrMtrcs_PregSplitEigRatio = nan; CheckErrMtrcs_PregSplit=false; end



%%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% do the partitioning
[ni,nj,~,~] = FindkPartition_Gander(A,SM_nmbsubdom_PwrOfTwo); rows_cumsum = cumsum(ni); cols_cumsum = cumsum(nj);
frst_rows = NaN(SM_nmbsubdom,1); last_rows = NaN(SM_nmbsubdom,1);  frst_cols = NaN(SM_nmbsubdom,1); last_cols = NaN(SM_nmbsubdom,1);
if strcmp(SM_type,'RAS'), mid_rows = NaN(SM_nmbsubdom,1); mid_cols = NaN(SM_nmbsubdom,1); end
for ind_subdom = 1:SM_nmbsubdom
    if ind_subdom == 1
        frst_rows(1) = 1; frst_cols(1) = 1;
    else
        frst_rows(ind_subdom) = rows_cumsum( (ind_subdom-1-1)*3 +1 ) + 1; frst_cols(ind_subdom) = cols_cumsum( (ind_subdom-1-1)*3 +1 ) + 1;
    end

    if strcmp(SM_type,'RAS') % if we do RAS -> get the midpoints of the subdomain overlaps
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
if debug == 2, plotPartition_Gander( A, ni, nj ); end % check for indices

if strcmp(SM_type,'RAS') % if we do RAS -> get the restriction indices
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






% %%% run Schwarz method
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_dd = u_init; res_prev = rhs-A*u_dd; 


if debug ~= 0 % track convergence
    u_exact = A\rhs; u_err = NaN(SM_nmbiter+1,N); 
    err_norm = NaN(SM_nmbiter+1,1); u_err(1,:) = u_exact - u_init; err_norm(1) = norm(u_err(1,:));
end

if plot_convergence % plot convergence
    %%% plot the approximate solution
    nmb_int_gridcols = length(x_mesh)-2; u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_exact,nmb_int_gridcols,nmb_int_gridcols)';
    mesh(x_mesh,x_mesh,u_plot); xlabel('x');ylabel('y'); view(angle1,angle2); title('Exact Solution'); pause(2);

    %%% plot the error
    nmb_int_gridcols = length(x_mesh)-2; u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_err(1,:),nmb_int_gridcols,nmb_int_gridcols)';
    mesh(x_mesh,x_mesh,u_plot); xlabel('x');ylabel('y'); view(angle1,angle2); title('Error'); pause(.1);
end


for ind_iter = 1:SM_nmbiter
    for ind_subdom = 1:SM_nmbsubdom
    
        %%% get my subdomain solutions  
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
                if strcmp('fp64',ChopOptions.format)
                    ExtraFactorForRHS = 100;
                elseif strcmp('fp32',ChopOptions.format)
                    ExtraFactorForRHS = 100;
                elseif strcmp('fp16',ChopOptions.format)
                    ExtraFactorForRHS = 100;
                elseif strcmp('bfloat16',ChopOptions.format)
                    ExtraFactorForRHS = 10;
                elseif strcmp('q43',ChopOptions.format)
                    ExtraFactorForRHS = 5;
                elseif strcmp('q52',ChopOptions.format)
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
                %uaux_check = curr_Lfact \ curr_rhs_chopped; norm(uaux-uaux_check)/norm(uaux_check)
                curr_u_bfrPerm_bfrScale = trisol_SprsMtrxDnseRhs(curr_Ufact,uaux, ChopOptions,true,curr_Ufact_NnzStruct);
                %u_check = curr_Ufact \ uaux; norm(curr_u_bfrPerm_bfrScale-u_check)/norm(u_check)

                % permute the solution vector back
                curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{9}; curr_u_Permed = double(curr_u_bfrPerm_bfrScale(curr_Qvec_inv));
                curr_u_Permed_Scaled = curr_u_Permed;
            end

            %%% Checking the theoretical requirements
            if ind_iter == 1
                if CheckErrMtrcs_norms % compute the norm || mathcalA_i^{-1}F_i ||
                    curr_ErrMtrx = MySubdomFactors_LP{ind_subdom}{10}; mathcalAi = MySubdomFactors_LP{ind_subdom}{12};
                    Ai_inv_ErrMtrx_Matvec = @(v,tflag) MatVecFun(v,tflag,mathcalAi,curr_ErrMtrx); 
                    ErrMtrcs_NormVals(ind_subdom) = svds(Ai_inv_ErrMtrx_Matvec,size(mathcalAi),1);
                    disp( "      For chop precision " + RoundingOptions{1}.format + " & i = " + ind_subdom + " we get norm(Ai_inv * Err_i) = " + ErrMtrcs_NormVals(ind_subdom))
                    if ErrMtrcs_NormVals(ind_subdom) < 1,  ErrMtrcs_NormCond(ind_subdom) = true; end
                end
                if CheckErrMtrcs_entrs % check the inequality "mathcalA_i^{-1} >= mathcalA_i^{-1} * F_i * mathcalA_i^{-1}" -- only for small enough sizes
                    curr_ErrMtrx = MySubdomFactors_LP{ind_subdom}{10}; mathcalAi = MySubdomFactors_LP{ind_subdom}{12};
                    Ainv = inv(mathcalAi); AinvErrAinv = Ainv * curr_ErrMtrx * Ainv;
                    indsOfTheViolation = Ainv <= AinvErrAinv; ErrMtrcs_EntrsViolation(ind_subdom) = norm( Ainv(indsOfTheViolation) - AinvErrAinv(indsOfTheViolation) );
                    disp( "      For chop precision " + RoundingOptions{1}.format + " & i = " + ind_subdom + " we get norm(Ai_inv - Ai_inv*Err_i*Ai_inv) = " + ErrMtrcs_EntrsViolation(ind_subdom))
                    if ErrMtrcs_EntrsViolation(ind_subdom) == 0,  ErrMtrcs_EntrCond(ind_subdom) = true; end
                end
                if CheckErrMtrcs_PregSplit % check the inequality "lambda_min(A) >= lambda_max(mu_niv*Di_inv*Err_i*Di_inv) -- only for small enough sizes
                    curr_ErrMtrx = MySubdomFactors_LP{ind_subdom}{10}; mathcalAi = MySubdomFactors_LP{ind_subdom}{12};
                    lambda_min = abs(eigs( mathcalAi, 1, 'smallestabs')); lambda_minInf_abs = abs(eigs( curr_ErrMtrx, 1, 'smallestreal')); ErrMtrcs_PregSplitEigRatio(ind_subdom) = lambda_min/lambda_minInf_abs;
                    disp( "      For chop precision " + RoundingOptions{1}.format + " & i = " + ind_subdom + " we get:" + newline + "        LambdaMin(mathcalAi) = " + lambda_min + newline + "        | Lambda_MinInf(Err_i) | = " + lambda_minInf_abs)
                    if lambda_min >= lambda_minInf_abs,  ErrMtrcs_PregSplitCond(ind_subdom) = true; end
                end
                % we do not check "F_i >= 0" as that is by construction
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

            %%% Checking the theoretical requirements
            if ind_iter == 1
                if CheckErrMtrcs_norms % compute the norm || mathcalA_i^{-1}F_i ||
                    curr_ErrMtrx = MySubdomFactors_LP{ind_subdom}{8}; mathcalAi = MySubdomFactors_LP{ind_subdom}{10};
                    Ai_inv_ErrMtrx_Matvec = @(v,tflag) MatVecFun(v,tflag,mathcalAi,curr_ErrMtrx); 
                    ErrMtrcs_NormVals(ind_subdom) = svds(Ai_inv_ErrMtrx_Matvec,size(mathcalAi),1);
                    disp( "      For " + Advanpix_NmbDigits + " digits & i = " + ind_subdom + " we get norm(Ai_inv * Err_i) = " + ErrMtrcs_NormVals(ind_subdom))
                    if ErrMtrcs_NormVals(ind_subdom) < 1,  ErrMtrcs_NormCond(ind_subdom) = true; end
                end
                if CheckErrMtrcs_entrs % check the inequality "mathcalA_i^{-1} >= mathcalA_i^{-1} * F_i * mathcalA_i^{-1}" -- only for small enough sizes
                    curr_ErrMtrx = MySubdomFactors_LP{ind_subdom}{8}; mathcalAi = MySubdomFactors_LP{ind_subdom}{10};
                    Ainv = inv(mathcalAi); AinvErrAinv = Ainv * curr_ErrMtrx * Ainv;
                    indsOfTheViolation = Ainv <= AinvErrAinv; ErrMtrcs_EntrsViolation(ind_subdom) = norm( Ainv(indsOfTheViolation) - AinvErrAinv(indsOfTheViolation) );
                    disp( "      For " + Advanpix_NmbDigits + " digits & i = " + ind_subdom + " we get norm(Ai_inv - Ai_inv*Err_i*Ai_inv) = " + ErrMtrcs_EntrsViolation(ind_subdom))
                    if ErrMtrcs_EntrsViolation(ind_subdom) <= 1e-20,  ErrMtrcs_EntrCond(ind_subdom) = true; end
                end
                if CheckErrMtrcs_PregSplit % check the inequality "lambda_min(A) >= lambda_max(mu_niv*Di_inv*Err_i*Di_inv) -- only for small enough sizes
                    curr_ErrMtrx = MySubdomFactors_LP{ind_subdom}{8}; mathcalAi = MySubdomFactors_LP{ind_subdom}{10};
                    lambda_min = abs(eigs( mathcalAi, 1, 'smallestabs')); lambda_minInf_abs = abs(eigs( curr_ErrMtrx, 1, 'smallestreal')); ErrMtrcs_PregSplitEigRatio(ind_subdom) = lambda_min/lambda_minInf_abs;
                    disp( "      For " + Advanpix_NmbDigits + " digits & i = " + ind_subdom + " we get:" + newline + "        LambdaMin(mathcalAi) = " + lambda_min + newline + "        | Lambda_MinInf(Err_i) | = " + lambda_minInf_abs)
                    if lambda_min >= lambda_minInf_abs,  ErrMtrcs_PregSplitCond(ind_subdom) = true; end
                end
                % we do not check "F_i >= 0" as that is by construction
            end
        end

        % %%% check
        % Ai = A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom));
        % curr_Qvec(curr_Qvec_inv) = 1:length(curr_Qvec_inv); %q_inv(q) = 1:length(q);
        % curr_u_check = Ai \ curr_rhs_bfrPerm_bfrScale; norm(curr_u_Permed_Scaled-curr_u_check)/norm(curr_u_check) %* 1/condest(Ai)
        % disp('ahoj')
        



        %%% the subdomain updating strategies, depending on Schwarz method type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ind_subdom == 1 %%% update my global solution for subdomain 1
            if strcmp(SM_type,'RAS') 
                u_dd = u_dd + [curr_u_Permed_Scaled(1:last_rows_RAS(ind_subdom)); zeros(N-last_rows_RAS(ind_subdom),1)];
            elseif strcmp(SM_type,'dAS')
                u_dd = u_dd + damping_theta * [curr_u_Permed_Scaled; zeros(N-last_rows(ind_subdom),1)];
            elseif strcmp(SM_type,'MS')
                u_dd = u_dd + [curr_u_Permed_Scaled; zeros(N-last_rows(ind_subdom),1)];
            end

        elseif ind_subdom == SM_nmbsubdom %%% update my global solution for subdomain N
            if strcmp(SM_type,'RAS') 
                frst_row_inside_subdom = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1;
                u_dd = u_dd + [zeros(last_rows_RAS(ind_subdom-1),1); curr_u_Permed_Scaled(frst_row_inside_subdom:end)]; 
            elseif strcmp(SM_type,'dAS')
                u_dd = u_dd + damping_theta * [zeros(frst_rows(ind_subdom)-1,1); curr_u_Permed_Scaled]; 
            elseif strcmp(SM_type,'MS')
                u_dd = u_dd + [zeros(frst_rows(ind_subdom)-1,1); curr_u_Permed_Scaled]; 
            end

        else %%% update my global solution for middle subdomains
            if strcmp(SM_type,'RAS') 
                frst_row_inside_subdom = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1; % nmb_of_rows_in_left_overlap = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) -> doesn't include the +1 to start after the overlap
                last_row_inside_subdom = last_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1; % nmb_of_rows_to_right_overlap = last_rows_RAS(ind_subdom) - frst_rows(ind_subdom) -> doesn't include the +1 to start after the overlap
                u_dd = u_dd + [zeros(last_rows_RAS(ind_subdom-1),1); curr_u_Permed_Scaled(frst_row_inside_subdom:last_row_inside_subdom); zeros(N-last_rows_RAS(ind_subdom),1)];
            elseif strcmp(SM_type,'dAS')
                u_dd = u_dd + damping_theta * [zeros(frst_rows(ind_subdom)-1,1); curr_u_Permed_Scaled; zeros(N-last_rows(ind_subdom),1)];
            elseif strcmp(SM_type,'MS')
                u_dd = u_dd + [zeros(frst_rows(ind_subdom)-1,1); curr_u_Permed_Scaled; zeros(N-last_rows(ind_subdom),1)];
            end
        end



        if strcmp(SM_type,'MS') % if multiplicative Schwarz -> update the residual after each subdomain solve
            res_prev = rhs-A*u_dd; 
        end
    end
    
    if strcmp(SM_type,'RAS') || strcmp(SM_type,'dAS') % if (restrcited) additive Schwarz -> update the residual after all subdomain solve
        res_prev = rhs-A*u_dd; 
    end
    u_seq(ind_iter,:) = u_dd; res_norm = norm(res_prev); 
    if res_norm/norm(rhs) < SM_relresacc, break; end

    if debug == 1 % track convergence 
        u_err(ind_iter+1,:) = u_exact - u_dd; err_norm(ind_iter+1) = norm(u_err(ind_iter+1,:)); 
        % if err_norm(ind_iter+1) < 1e-15 && err_norm(ind_iter) < 1e-15
        %     break
        % end
    end
    if plot_convergence % plot convergence
        % %%% plot the approximate solution
        % u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_dd,nmb_int_gridcols,nmb_int_gridcols)'; mesh(x_mesh,x_mesh,u_plot); view(angle1,angle2); pause(waittime);
        %%% plot the error
        u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_err(ind_iter+1,:),nmb_int_gridcols,nmb_int_gridcols)'; mesh(x_mesh,x_mesh,u_plot); xlabel('x');ylabel('y'); view(angle1,angle2); title('Error'); pause(waittime);
    end
end


MyOutput = {u_seq};
if debug == 1 % track convergence 
    if plot_convergence, semilogy([0,1:SM_nmbiter], err_norm/err_norm(1), 'ro-'); pause(1); end
    MyOutput = {u_seq,u_exact,err_norm , ErrMtrcs_NormCond,ErrMtrcs_NormVals , ErrMtrcs_EntrCond,ErrMtrcs_EntrsViolation , ErrMtrcs_PregSplitCond,ErrMtrcs_PregSplitEigRatio};
end






end