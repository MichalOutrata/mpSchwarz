clear; clc; close('all');

%%% choose set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Schwarz set-up
SM_methods = {"dAS", "RAS", "MS"}; % which methods we want to run out of "dAS", "RAS", "MS" 
SM_nmbsubdoms_PwrOfTwo = 1; % number of subdomains for the Schwarz method
SM_nmbiter = 42; % number of RAS iterations
SM_relresacc = 42.42; % number of RAS iterations
dampingTheta = 1/3; % generally one should use 1/(2^SM_nmbsubdoms_PwrOfTwo+1) but since we have sausage-like domain decomposition, we can do with alternating coloring no matter how many subdomains
SM_RandomRHS = true;
SchwarzMethodStuff = {SM_methods, 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc, dampingTheta};

%%% RoundingStuff set-up
RM_Software = "advanpix";   % "advanpix" or "chop"
RM_Type = "Mmtrx";          % "Mmtrx" or "StieltjessDiagDom"
if strcmp(RM_Software,"advanpix")
    RM_Precision_list = 1:16; % number of digits advanpix uses for the subdomain solves
elseif strcmp(RM_Software,"chop")
    RM_Precision_list = {'q52','q43','bfloat16','fp16','fp32','fp64'}; % precision format of chop for the subdomain solves
    ChopOptions.explim = 1; % do we consider the "range limitation". If "0" then bound on the maximal exponent is ignored inside chop   
    RM_ChopOptions = ChopOptions;
end
RM_ScalingToRange = true;           % whether or not we rescaled the matrix to fit to the considered low-precision range of representation 
RM_RescaleType = "diag2side_symm";       % "diag2side" or "diag2side_symm"
RM_ScalingToRange_Fraction = 16;    % we rescale the matrix and RHS to fit to "x_max / RM_ScalingToRange_Fraction"
RM_ScalingOptions = {RM_ScalingToRange,RM_RescaleType,RM_ScalingToRange_Fraction};


%%% choose GMRES set-up
GMRES_PrecType = 'L'; % type of preconditioner - none/left/right
GMRES_nmbiter = 100; % number of GMRES iterations
GMRES_relresacc = 1e-12; % GMRES accuracy
Prec = @mpSchwarzMethodsPreconds;

debug = 0;

%%% which problems we will run (see below)
list_of_nmb_int_gridcols = 50:20:330; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
if strcmp(RM_Type,"Mmtrx") && strcmp(RM_RescaleType,"diag2side") 
    ProblemChoices_list = [1,2,3]; SpecialStrngForNamingOutput = "NonSymmPrblms_N" + list_of_nmb_int_gridcols(end);% number of digits advanpix uses for the subdomain solves
elseif strcmp(RM_Type,"Mmtrx") && strcmp(RM_RescaleType,"diag2side_symm") 
    ProblemChoices_list = [4,5,6]; SpecialStrngForNamingOutput = "SymmPrblms_N" + list_of_nmb_int_gridcols(end); % number of digits advanpix uses for the subdomain solves
elseif strcmp(RM_Type,"StieltjessDiagDom") 
    ProblemChoices_list = [4,5,6]; SpecialStrngForNamingOutput = "SymmPrblms_DiagDomm_N" + list_of_nmb_int_gridcols(end);  % number of digits advanpix uses for the subdomain solves
end


%%% initializing the arrays for saving convergence
GMRESnoprec_ConvCurves = zeros( length(ProblemChoices_list), length(list_of_nmb_int_gridcols), GMRES_nmbiter+1 );
GMRESnoprec_NmbItToConv = zeros( length(ProblemChoices_list), length(list_of_nmb_int_gridcols) );
GMRESprec_ConvCurves = zeros( length(ProblemChoices_list), length(list_of_nmb_int_gridcols), length(RM_Precision_list), length(SM_methods), GMRES_nmbiter+1 );
GMRESprec_NmbItToConv = zeros( length(ProblemChoices_list), length(list_of_nmb_int_gridcols), length(RM_Precision_list), length(SM_methods) );






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
        % rhs = zeros(size(A,1),1); rhs(1:nmb_int_gridcols:end)=-(-1)/h^2; % right-hand side, one minus for shifting to the RHS, one minus for having negative laplacian
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
        rhs = rand(size(A,1),1)*1/h^2; GMRES_initguess = zeros(size(rhs)); 
    else
        GMRES_initguess = zeros(size(A,1),1); rhs = ones(size(A,1),1); % rhs = ones(size(A,1),1)*1/h^2;
    end   


    %%% run non-preconditioned GMRES
    PrecInfo = {'noprec'}; disp("   No preconditioner")
    [AddOutputs,x, res_norms] = MyGMRES_PrecApply(A, rhs, PrecInfo, Prec, GMRES_initguess, GMRES_nmbiter, GMRES_relresacc );
    if AddOutputs{4} >= 1 
        conv_curves = res_norms/res_norms(1); conv_curves_prettier = nan(GMRES_nmbiter+1,1); conv_curves_prettier(1:length(conv_curves)) = conv_curves;
        for ind = 1:GMRES_nmbiter
            if ind <= length(conv_curves) && conv_curves(ind) < GMRES_relresacc, conv_curves_prettier(ind:end) = nan; break; end
        end
    end
    GMRESnoprec_ConvCurves(ind_Prblm, ind_meshsize, :) = conv_curves_prettier;
    GMRESnoprec_NmbItToConv(ind_Prblm, ind_meshsize) = AddOutputs{4};



    for ind_whichPrecision = 1:length(RM_Precision_list)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(RM_Software,"advanpix")
            disp("   " + "advanpix d_s = " + RM_Precision_list(ind_whichPrecision))
            RM_PrecisionChoice = RM_Precision_list(ind_whichPrecision);
        elseif strcmp(RM_Software,"chop")
                disp("   " + " chop format = " + RM_Precision_list{ind_whichPrecision})
            RM_ChopOptions.format = RM_Precision_list{ind_whichPrecision}; RM_PrecisionChoice = RM_ChopOptions;
        end
        RM_RoundingOptions = {RM_PrecisionChoice,RM_ScalingOptions};
        RoundingStuff = {RM_Software, RM_Type, RM_RoundingOptions};
        [MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethodStuff,RoundingStuff, debug);

        


        for ind_whichMethod = 1:length(SM_methods)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp("     - " +  SM_methods{ind_whichMethod})

            MyPrecondPackage{1}{1} = SM_methods{ind_whichMethod};
            PrecInfo = {GMRES_PrecType,'Handle',{A,MyPrecondPackage}};
            [AddOutputs,x, res_norms] = MyGMRES_PrecApply(A, rhs, PrecInfo, Prec, GMRES_initguess, GMRES_nmbiter, GMRES_relresacc );

            if AddOutputs{8} >= 1 
                conv_curves = res_norms/res_norms(1); conv_curves_prettier = nan(GMRES_nmbiter+1,1); conv_curves_prettier(1:length(conv_curves)) = conv_curves;
                for ind = 1:GMRES_nmbiter
                    if ind <= length(conv_curves) && conv_curves(ind) < GMRES_relresacc, conv_curves_prettier(ind:end) = nan; break; end
                end
            end
            GMRESprec_ConvCurves(ind_Prblm, ind_meshsize, ind_whichPrecision, ind_whichMethod,:) = conv_curves_prettier;
            GMRESprec_NmbItToConv(ind_Prblm, ind_meshsize, ind_whichPrecision, ind_whichMethod) = AddOutputs{8};
        end
    end
end
end



%%% Packaging
Info_SchwarzMethodStuff = {'which Schwarz methods we used',SM_methods,'# subdoms = 2^(...)',SM_nmbsubdoms_PwrOfTwo,'theta for dAS', dampingTheta, 'use random RHS', SM_RandomRHS}; 
Info_GMRESStuff = {'GMRES preconditioner - none/left/right',GMRES_PrecType,'# GMRES it,',GMRES_nmbiter,'GMRES rel. res. acc', GMRES_relresacc}; 
Info_RoundingStuff = {'used rounding software',RM_Software, 'used rounding approach - Mmtrx vs Stieltjes',RM_Type, 'which precisions we used',RM_Precision_list, 'what scaling set-up we used',RM_ScalingOptions };
Info_PrblmAndDiscr = {'problem indices', ProblemChoices_list, '# interior grid points on the side of the square',list_of_nmb_int_gridcols};
Output_Convergence = {'Preconditioned GMRES ConvCurvs', GMRESprec_ConvCurves, '# iterations for preconditioned GMRES to converge', GMRESprec_NmbItToConv, 'NON-Preconditioned GMRES ConvCurvs', GMRESnoprec_ConvCurves, '# iterations for NON-preconditioned GMRES to converge', GMRESnoprec_NmbItToConv };

MyData = {'specs used for Schwarz method',Info_SchwarzMethodStuff, ... 
    'specs used for GMRES',Info_GMRESStuff, ... 
    'specs used for rounding/multiprec',Info_RoundingStuff, ... 
    'specs for discretization and problem generation', Info_PrblmAndDiscr, ...
    'outputs for convergence curves and convergence factors', Output_Convergence};

s_SaveString = "SavedData_mpSMprec_" + SpecialStrngForNamingOutput + "_" + RM_Software + "_Scaling" + RM_ScalingToRange + "_RandomRHS" + SM_RandomRHS + ".mat";
save(s_SaveString,'MyData');  


%%% plot the ConvFact as a function of number of digits
ind_Prblm = 2; % choose
ind_meshsize = 2; % choose
ind_whichMethod = 2; % choose

ClsstIntSqrt = floor(sqrt(length(RM_Precision_list))); nmbColsTiles = floor( length(RM_Precision_list) / ClsstIntSqrt );
if length(RM_Precision_list) <= 3, nmbColsTiles = length(RM_Precision_list); nmbRowsTiles = 1; else
    if floor( length(RM_Precision_list) / nmbColsTiles ) == length(RM_Precision_list) / nmbColsTiles, nmbRowsTiles = floor( length(RM_Precision_list) / nmbColsTiles ); 
    else, nmbRowsTiles = floor( length(RM_Precision_list) / nmbColsTiles ) + 1; end
end
figure(96); t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');
for ind_nmbdig = 1:length(RM_Precision_list)
    ax(ind_nmbdig) = nexttile();
    PlotData = reshape( GMRESprec_ConvCurves(ind_Prblm, ind_meshsize, ind_nmbdig, ind_whichMethod,:), GMRES_nmbiter+1,1);
    semilogy(PlotData,'ko-');
end
title(t, append('# grid points on one side = ',num2str(list_of_nmb_int_gridcols(ind_meshsize))));

figure(97)
plot(RM_Precision_list, reshape( GMRESprec_NmbItToConv(ind_Prblm, ind_meshsize, :, ind_whichMethod), size(RM_Precision_list) ),'ro-'); xlabel('# digits');ylabel('# GMRES it to conv'); title( append('# grid points on one side = ',num2str(list_of_nmb_int_gridcols(ind_meshsize))));