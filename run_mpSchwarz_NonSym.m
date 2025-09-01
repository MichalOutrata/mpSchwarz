clear; clc; close('all');

%%% choose set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Schwarz set-up
SM_methods = {"dAS", "RAS", "MS"}; % which methods we want to run out of "dAS", "RAS", "MS" 
SM_nmbsubdoms_PwrOfTwo = 1; % number of subdomains for the Schwarz method
SM_nmbiter = 61; % number of RAS iterations
SM_relresacc = 1e-12; % number of RAS iterations
dampingTheta = 1/3; % generally one should use 1/(2^SM_nmbsubdoms_PwrOfTwo+1) but since we have sausage-like domain decomposition, we can do with alternating coloring no matter how many subdomains
SM_RandomRHS = true;

%%% RoundingStuff set-up
RM_Software = "chop";   % "advanpix" or "chop"
RM_Type = "Mmtrx";          % "Mmtrx" or "StieltjessDiagDom"
if strcmp(RM_Software,"advanpix")
    RM_Precision_list = 1:16; % number of digits advanpix uses for the subdomain solves
elseif strcmp(RM_Software,"chop")
    RM_Precision_list = {'q52','q43','bfloat16','fp16','fp32','fp64'}; %{'q43','bfloat16','fp16'}; % precision format of chop for the subdomain solves
    ChopOptions.explim = 1; % do we consider the "range limitation". If "0" then bound on the maximal exponent is ignored inside chop   
    RM_ChopOptions = ChopOptions;
end
RM_ScalingToRange = true;           % whether or not we rescaled the matrix to fit to the considered low-precision range of representation 
RM_RescaleType = "diag2side";       % "diag2side" or "diag2side_symm"
RM_ScalingToRange_Fraction = 16;    % we rescale the matrix and RHS to fit to "x_max / RM_ScalingToRange_Fraction"
RM_ScalingOptions = {RM_ScalingToRange,RM_RescaleType,RM_ScalingToRange_Fraction};

%%% ErrorCheckStuff set-up
CheckErrMtrcs_norms = true; % whether or not we check the condition ||mathcalAi_inv * F_i || < 1  
CheckErrMtrcs_entrs = true; % whether or not we check the condition mathcalAi_inv >= mathcalAi_inv * F_i * mathcalAi_inv   
CheckErrMtrcs_PregSplit = false; % whether or not we check the condition mathcalAi_inv >= mathcalAi_inv * F_i * mathcalAi_inv (considered only for "StieltjessDiagDom")    
if strcmp(RM_Type,"Mmtrx") && strcmp(RM_RescaleType,"diag2side"), CheckErrMtrcs_PregSplit = false; end

%%% which problems we will run (see below)
list_of_nmb_int_gridcols = 50:5:55; %50:20:330; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
if strcmp(RM_Type,"Mmtrx") && strcmp(RM_RescaleType,"diag2side") 
    ProblemChoices_list = [1,2,3]; SpecialStrngForNamingOutput = "NonSymmPrblms_N" + list_of_nmb_int_gridcols(end);% number of digits advanpix uses for the subdomain solves
elseif strcmp(RM_Type,"Mmtrx") && strcmp(RM_RescaleType,"diag2side_symm") 
    ProblemChoices_list = [4,5,6]; SpecialStrngForNamingOutput = "SymmPrblms_N" + list_of_nmb_int_gridcols(end); % number of digits advanpix uses for the subdomain solves
elseif strcmp(RM_Type,"StieltjessDiagDom") 
    ProblemChoices_list = [4,5,6]; SpecialStrngForNamingOutput = "SymmPrblms_DiagDomm_N" + list_of_nmb_int_gridcols(end);  % number of digits advanpix uses for the subdomain solves
end

%%% processing results set-up
DoPlotting = false; 
SavePlotData = true;

%%% plotting
indsPrblms_PlotErr = 1:length(ProblemChoices_list); % we save data for plotting for all considered problems
indsMesh_PlotErr = 1; % we save data for plotting ONLY for the smallest resolution
indsIter_PlotErr = [1,2,3]; % we save iterations 10, 20 and 60 for plotting
if strcmp(RM_Software,"advanpix")
    indsWhichPrecs_PlotErr = [2,3,4]; counter_WhichPrecPlotErr = 1;
elseif strcmp(RM_Software,"chop")
    indsWhichPrecs_PlotErr = [2,3,4]; counter_WhichPrecPlotErr = 1;
end
PlotData_ErrPlot = cell(length(ProblemChoices_list),length(SM_methods),length(indsWhichPrecs_PlotErr),length(indsIter_PlotErr));
u_ExactSols = cell(length(ProblemChoices_list));

%%% initializing the arrays for saving convergence & theory results 
if CheckErrMtrcs_norms, ErrMtrcsNrms_CondStsfd = false(length(ProblemChoices_list),length(list_of_nmb_int_gridcols),length(RM_Precision_list)); else, ErrMtrcsNrms_CondStsfd = nan; end
if CheckErrMtrcs_entrs, ErrMtrcsEntrs_CondStsfd = false(length(ProblemChoices_list),length(list_of_nmb_int_gridcols),length(RM_Precision_list)); else, ErrMtrcsEntrs_CondStsfd = nan; end
if CheckErrMtrcs_PregSplit, ErrMtrcsPregSplits_CondStsfd = false(length(ProblemChoices_list),length(list_of_nmb_int_gridcols),length(RM_Precision_list)); else, ErrMtrcsPregSplits_CondStsfd = nan; end
ConvCurves = zeros( length(ProblemChoices_list), length(list_of_nmb_int_gridcols), length(SM_methods), length(RM_Precision_list), SM_nmbiter+1 );
ConvFactApprox = zeros( length(ProblemChoices_list), length(list_of_nmb_int_gridcols), length(SM_methods), length(RM_Precision_list) );






for ind_Prblm = 1:length(ProblemChoices_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
ProblemChoice = ProblemChoices_list(ind_Prblm);

for ind_meshsize = 1:length(list_of_nmb_int_gridcols)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    nmb_int_gridcols = list_of_nmb_int_gridcols(ind_meshsize); h=1/(nmb_int_gridcols+1); % mesh size
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp("Problem " + ProblemChoice + " of size = " + nmb_int_gridcols^2)


    if ProblemChoice == 0 %%% do negative laplacian
        G=numgrid('S',nmb_int_gridcols+2);
        eta = @(x,y) zeros(size(x)); a = @(x,y) ones(size(x)); b1 = @(x,y) zeros(size(x)); b2 = b1;
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,a,b1,b2});
        rhs = zeros(size(A,1),1); rhs(1:nmb_int_gridcols:end)=-(-1); rhs = 1/h^2*rhs; % right-hand side, one minus for shifting to the RHS, one minus for having negative laplacian
        u_plot = zeros(nmb_int_gridcols+2); %u_plot(2:nmb_int_gridcols+1,1) = 1; 
        angle1 = 52.5; angle2 = 30; x_mesh = 0:h:1;

    elseif ProblemChoice == 1 %%% do non-symmetric AdvecDiff based on SzyldGander Fig 2.1
        G=numgrid('S',nmb_int_gridcols+2);
        eta = @(x,y) x.^2.*cos(x+y).^2'; a = @(x,y) 20*(x+y).^2.*exp(x-y); b1 =  @(x,y) (y-0.5); b2 =  @(x,y) -(x-0.5);
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,a,b1,b2});
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; x_mesh = 0:h:1; 
    
    elseif ProblemChoice == 2 %%% do non-symmetric AdvecDiff based SzyldFrommer eqn (15) and p.660 -> changed the advection strength to ~40 so that we have an Mmtrx
        G=numgrid('S',nmb_int_gridcols+2); AdvectStrngth = 100;
        eta = @(x,y) 0.*x+0.*y; alpha = @(x,y) 1 + 0.*x+0.*y; beta = @(x,y) 1 + 0.*x+0.*y; nu = @(x,y) AdvectStrngth.*(1.*x.*(x-1).*(1-2.*y)); mu = @(x,y) AdvectStrngth.*(-1.*y.*(y-1).*(1-2.*x));
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,alpha,beta,nu,mu}); %A_szyld = A*h^2+speye(size(A,1)); disp(isMmtrx(A)); disp(isMmtrx(A_szyld)); disp(issymmetric(A))
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; x_mesh = 0:h:1;

    elseif ProblemChoice == 3 %%% based on Problem 2, only we added a region of much higher diffusion to artificially bump up condition number
        G=numgrid('S',nmb_int_gridcols+2);
        AdvectStrngth = 100;
        clearvars alpha eta b1 b2
        addpath(genpath('./CoeffFuncs_ReacAdvecDiff'));
        eta_handle = @(x,y) eta(x,y); a_handle = @(x,y) alpha(x,y); b1_handle = @(x,y) b1(x,y,AdvectStrngth); b2_handle = @(x,y) b2(x,y,AdvectStrngth);
        A=ReactAdvDiff_Sqr_FD('SzyldGander_SclrFeval',G,{eta_handle,a_handle,b1_handle,b2_handle});
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; x_mesh = 0:h:1; 


    elseif ProblemChoice == 4 %%% do symmetric AdvecDiff based on SzyldGander eqn (15,18) p.648
        G=numgrid('S',nmb_int_gridcols+2);
        eta = @(x,y) x.^2.*cos(x+y).^2'; alpha =  @(x,y) (x+y).^2.*exp(x-y); nu = @(x,y) 0*(1.*x+1.*y); mu = @(x,y) 0*(1.*x+1.*y);
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,alpha,nu,mu});
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; x_mesh = 0:h:1;
    
    elseif ProblemChoice == 5 %%% do symmetric AdvecDiff based on SzyldFrommer eqn (15,18) p.648
        G=numgrid('S',nmb_int_gridcols+2); 
        eta = @(x,y) 500.*x+1.*y; alpha = @(x,y) 1+9.*(x+y); beta = @(x,y) 1+9.*(x+y); nu = @(x,y) 0.*x+0.*y; mu = @(x,y) 0.*x+0.*y;
        A=ReactAdvDiff_Sqr_FD('SzyldFrommer',G,{eta,alpha,beta,nu,mu});
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; x_mesh = 0:h:1;

    elseif ProblemChoice == 6 %%% based on Problem 4, only we added a region of much higher diffusion to artificially bump up condition number
        G=numgrid('S',nmb_int_gridcols+2);
        AdvectStrngth = 0;
        clearvars alpha eta b1 b2
        addpath(genpath('./CoeffFuncs_ReacAdvecDiff'));
        eta_handle = @(x,y) eta(x,y); a_handle = @(x,y) alpha(x,y); b1_handle = @(x,y) b1(x,y,AdvectStrngth); b2_handle = @(x,y) b2(x,y,AdvectStrngth);
        A=ReactAdvDiff_Sqr_FD('SzyldGander_SclrFeval',G,{eta_handle,a_handle,b1_handle,b2_handle});
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; x_mesh = 0:h:1;

    % %%% optional more advanced rhs building
    % [m,n] = size(G); h=1/(n-1); [X,Y]=meshgrid(0:h:1,0:h:(m-1)*h); IntriorInds = find(G);
    % rhsfun = @(x,y) x.*y.*(1-x).*(1-y); rhsfun_discr = feval(rhsfun,X,Y); rhs = A*rhsfun_discr(IntriorInds);
    end

    if SM_RandomRHS
        rhs = rand(size(A,1),1); 
        u_init = -rand(size(A,1),1);
    else
        rhs = ones(size(A,1),1); % rhs = ones(size(A,1),1)*1/h^2;
        u_init = zeros(size(A,1),1);
    end
    if ind_meshsize == indsMesh_PlotErr, u_ExactSols{ind_Prblm} = A\rhs; end
    

    debug = 1; if DoPlotting, MyPlot = {x_mesh,u_plot,0.05,angle1,angle2}; else, MyPlot = {}; end % MyPlot = {x,u_plot,waitime,angle1,angle2};

    for ind_whichMethod = 1:length(SM_methods)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SM_method = SM_methods{ind_whichMethod};
        SchwarzMethod = {SM_method, 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc, dampingTheta}; 
        
        if ind_whichMethod == 1
            % we only need to check the conditions for the subdomain matrices, ie, the conditions are emthod-independent, hence we only check for the frst method 
            if CheckErrMtrcs_norms, CheckErrMtrcs_norms_StillCheck = true; else, CheckErrMtrcs_norms_StillCheck = false;  end
            if CheckErrMtrcs_entrs, CheckErrMtrcs_entrs_StillCheck = true; else, CheckErrMtrcs_entrs_StillCheck = false;  end
            if CheckErrMtrcs_PregSplit, CheckErrMtrcs_PregSplit_StillCheck = true; else, CheckErrMtrcs_PregSplit_StillCheck = false;  end
        else
            CheckErrMtrcs_norms_StillCheck = false; CheckErrMtrcs_entrs_StillCheck = false; CheckErrMtrcs_PregSplit_StillCheck = false;
        end 
        
        
        for ind_whichPrecision = 1:length(RM_Precision_list)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('   -------------------')
            if strcmp(RM_Software,"advanpix")
                disp("   " +  SM_method + " with advanpix d_s = " + RM_Precision_list(ind_whichPrecision))
                RM_PrecisionChoice = RM_Precision_list(ind_whichPrecision);
            elseif strcmp(RM_Software,"chop")
                disp("   " +  SM_method + " with chop format = " + RM_Precision_list{ind_whichPrecision})
                RM_ChopOptions.format = RM_Precision_list{ind_whichPrecision};
                RM_PrecisionChoice = RM_ChopOptions;
            end
            RM_RoundingOptions = {RM_PrecisionChoice,RM_ScalingOptions};
            RoundingStuff = {RM_Software, RM_Type, RM_RoundingOptions};
            ErrorCheckStuff = {CheckErrMtrcs_norms_StillCheck,CheckErrMtrcs_entrs_StillCheck,CheckErrMtrcs_PregSplit_StillCheck};
            [SMoutput] = mpSchwarzMethods(A,rhs,u_init, SchwarzMethod, RoundingStuff, ErrorCheckStuff, MyPlot, debug);
    
            %%% get the approximate convergence factors
            % the "SMoutput" is ordered as: SMoutput = {u_seq,u_exact,err_norm, ... 
            %     ErrMtrcs_NormCond,ErrMtrcs_NormVals , ErrMtrcs_EntrCond,ErrMtrcs_EntrsViolation , ErrMtrcs_PregSplitCond,ErrMtrcs_PregSplitEigRatio};
            error_nrms = SMoutput{3}; error_nrms = error_nrms / error_nrms(1); ConvCurves(ind_Prblm,ind_meshsize,ind_whichMethod,ind_whichPrecision,:) = error_nrms; 
            ind_aux = 1;
            while ind_aux < length(error_nrms) && error_nrms(ind_aux+1) > 1e-15, ind_aux = ind_aux + 1; end; ind_EndOfConv = ind_aux;
            ConsecErrsRatio = error_nrms(2:ind_EndOfConv)./error_nrms(1:ind_EndOfConv-1);
            ConvFactApproxSequence = sqrt( ConsecErrsRatio(1:end-1) .* ConsecErrsRatio(2:end) ); % correct smoothing of the ConvFactr due to having 2 subdomains nto perfectly symmetrical in space nor coeffs
            % disp( num2str( norm(ConvFactApproxSequence-ConvFactApproxSequence1)/norm(ConvFactApproxSequence) ))
            if ConvFactApproxSequence(end) > 1 
                ConvFactApprox(ind_Prblm,ind_meshsize,ind_whichMethod,ind_whichPrecision) = nan; ConvCurves(ind_Prblm,ind_meshsize,ind_whichMethod,ind_whichPrecision,:) = nan(size(error_nrms));
            else
                ConvFactApprox(ind_Prblm,ind_meshsize,ind_whichMethod,ind_whichPrecision) = ConvFactApproxSequence(end);
            end
 
            %%% get the errors when we want them
            % the "SMoutput" is ordered as: SMoutput = {u_seq,u_exact,err_norm, ... 
            %     ErrMtrcs_NormCond,ErrMtrcs_NormVals , ErrMtrcs_EntrCond,ErrMtrcs_EntrsViolation , ErrMtrcs_PregSplitCond,ErrMtrcs_PregSplitEigRatio};
            if ind_meshsize == indsMesh_PlotErr
                if indsWhichPrecs_PlotErr(counter_WhichPrecPlotErr) == ind_whichPrecision

                    if indsIter_PlotErr(1) == 0
                        u_sol = u_init; u_PlotErr = reshape( SMoutput{2}-u_sol, [length(rhs),1] );
                        if SavePlotData, PlotData_ErrPlot{ind_Prblm, ind_whichMethod, counter_WhichPrecPlotErr, 1} = u_PlotErr;  end
                        for ind = 2:length(indsIter_PlotErr)
                            u_solseq = SMoutput{1}; u_sol = u_solseq(indsIter_PlotErr(ind),:); 
                            u_PlotErr = reshape( SMoutput{2}'-u_sol, [length(rhs),1] );
                            if SavePlotData, PlotData_ErrPlot{ind_Prblm, ind_whichMethod, counter_WhichPrecPlotErr, ind} = u_PlotErr;  end
                        end
                        if counter_WhichPrecPlotErr < length(indsWhichPrecs_PlotErr), counter_WhichPrecPlotErr=counter_WhichPrecPlotErr+1; else, counter_WhichPrecPlotErr=1; end

                    else
                        for ind = 1:length(indsIter_PlotErr)
                            u_solseq = SMoutput{1}; u_sol = u_solseq(indsIter_PlotErr(ind),:); 
                            u_PlotErr = reshape( SMoutput{2}'-u_sol, [length(rhs),1] );
                            if SavePlotData, PlotData_ErrPlot{ind_Prblm, ind_whichMethod, counter_WhichPrecPlotErr, ind} = u_PlotErr;  end
                        end
                        if counter_WhichPrecPlotErr < length(indsWhichPrecs_PlotErr), counter_WhichPrecPlotErr=counter_WhichPrecPlotErr+1; else, counter_WhichPrecPlotErr=1; end

                    end

                end
            end
    
            %%% check the theory if we want to
            % the "SMoutput" is ordered as: SMoutput = {u_seq,u_exact,err_norm, ... 
            %     ErrMtrcs_NormCond,ErrMtrcs_NormVals , ErrMtrcs_EntrCond,ErrMtrcs_EntrsViolation , ErrMtrcs_PregSplitCond,ErrMtrcs_PregSplitEigRatio};
            if strcmp(RM_Software,"advanpix")
            % for advanpix, the precisions are ordered nicely, so
            % for runtime reasons, we keep checking the convergence conditions only until theyre satisfied for a certain precision for all subdomains 
                if CheckErrMtrcs_norms_StillCheck 
                    if SMoutput{4}, ErrMtrcsNrms_CondStsfd(ind_Prblm,ind_meshsize,ind_whichPrecision:end)=true; CheckErrMtrcs_norms_StillCheck=false; end
                end
                if CheckErrMtrcs_entrs_StillCheck
                    if SMoutput{6}, ErrMtrcsEntrs_CondStsfd(ind_Prblm,ind_meshsize,ind_whichPrecision:end)=true; CheckErrMtrcs_entrs_StillCheck=false; end
                end
                if CheckErrMtrcs_PregSplit_StillCheck
                    if SMoutput{8}, ErrMtrcsPregSplits_CondStsfd(ind_Prblm,ind_meshsize,ind_whichPrecision:end)=true; CheckErrMtrcs_PregSplit_StillCheck=false; end
                end
            elseif strcmp(RM_Software,"chop")
            % for chop, some consecutive precisions are comparable (q43&q52 or h&bfloat), so
            % we run everything till 'h' (assuming, 'h' is after 'bfloat16') and then keep checking only until theyre satisfied for a certain precision for all subdomains 
                if CheckErrMtrcs_norms_StillCheck 
                    if SMoutput{4}
                        ErrMtrcsNrms_CondStsfd(ind_Prblm,ind_meshsize,ind_whichPrecision)=true; 
                        if any(strcmp({'fp16','fp32','fp64'},RM_Precision_list{ind_whichPrecision}))
                            CheckErrMtrcs_norms_StillCheck=false; ErrMtrcsNrms_CondStsfd(ind_Prblm,ind_meshsize,ind_whichPrecision:end)=true;
                        end
                    end
                end
                if CheckErrMtrcs_entrs_StillCheck
                    if SMoutput{6}
                        ErrMtrcsEntrs_CondStsfd(ind_Prblm,ind_meshsize,ind_whichPrecision)=true; 
                        if any(strcmp({'fp16','fp32','fp64'},RM_Precision_list{ind_whichPrecision}))
                            CheckErrMtrcs_entrs_StillCheck=false; ErrMtrcsEntrs_CondStsfd(ind_Prblm,ind_meshsize,ind_whichPrecision:end)=true;
                        end                       
                    end
                end
                if CheckErrMtrcs_PregSplit_StillCheck
                    if SMoutput{8}
                        ErrMtrcsPregSplits_CondStsfd(ind_Prblm,ind_meshsize,ind_whichPrecision)=true; 
                        if any(strcmp({'fp16','fp32','fp64'},RM_Precision_list{ind_whichPrecision}))
                            CheckErrMtrcs_PregSplit_StillCheck=false; ErrMtrcsPregSplits_CondStsfd(ind_Prblm,ind_meshsize,ind_whichPrecision:end)=true;
                        end 
                    end
                end
            end


        end
    

    end
end
end


%%% Packaging
Info_SchwarzMethodStuff = {'which Schwarz methods we used',SM_methods,'# subdoms = 2^(...)',SM_nmbsubdoms_PwrOfTwo,'SM #iter', SM_nmbiter,'SM relres accuracy stop', SM_relresacc, 'theta for dAS', dampingTheta, 'use random RHS & initial guess', SM_RandomRHS}; 
Info_RoundingStuff = {'used rounding software',RM_Software, 'used rounding approach - Mmtrx vs Stieltjes',RM_Type, 'which precisions we used',RM_Precision_list, 'what scaling set-up we used',RM_ScalingOptions };
Info_ErrorCheckStuff = {'check cond. on norm',CheckErrMtrcs_norms, 'check cond. on entries of inv',CheckErrMtrcs_entrs, 'check cond. on P-regular splitting (eigvals)',CheckErrMtrcs_PregSplit};
Info_PrblmAndDiscr = {'problem indices', ProblemChoices_list, '# interior grid points on the side of the square',list_of_nmb_int_gridcols};
Output_Plotting = {'exact solution (only for the largest problem)', u_ExactSols, 'Error plots (only for the largest problem)', PlotData_ErrPlot, ...
    'we saved all problem indicies and only the coarsest (=first) mesh resolution','indices of accuracies/digits for which we plot error', indsWhichPrecs_PlotErr, 'SM iteration indices for which we plot error', indsIter_PlotErr,};
Output_Convergence = {'ConvCurvs for Schwarz methods', ConvCurves, 'approx ConvFacts for Schwarz methods', ConvFactApprox};
Output_TheoryConds = {'condition ||F_imathcalAiInv|| < 1',ErrMtrcsNrms_CondStsfd, 'condition mathcalAiInv >= mathcalAiInvF_imathcalAiInv',ErrMtrcsEntrs_CondStsfd, 'condition lambdamin(A_i) >= lambdaMinInf(Fi)',ErrMtrcsPregSplits_CondStsfd};

MyData = {'specs used for Schwarz method',Info_SchwarzMethodStuff, ... 
    'specs used for rounding/multiprec',Info_RoundingStuff, ... 
    'specs used for theory checks',Info_ErrorCheckStuff, ...
    'specs for discretization and problem generation', Info_PrblmAndDiscr, ...
    'outputs for plotting solutions and errors', Output_Plotting, ...
    'outputs for convergence curves and convergence factors', Output_Convergence, ...
    'outputs for checking theoretical conditions being satisfied', Output_TheoryConds};

s_SaveString = "SavedData_mpSM_" + SpecialStrngForNamingOutput + "_" + RM_Software + "_Scaling" + RM_ScalingToRange + "_RandomRHS" + SM_RandomRHS + ".mat";
save(s_SaveString,'MyData');  



