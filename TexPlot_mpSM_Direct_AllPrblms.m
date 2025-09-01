% Info_SchwarzMethodStuff = {'which Schwarz methods we used',SM_methods,'# subdoms = 2^(...)',SM_nmbsubdoms_PwrOfTwo,'SM #iter', SM_nmbiter,'SM relres accuracy stop', SM_relresacc, 'theta for dAS', dampingTheta}; 
% Info_RoundingStuff = {'used rounding software',RM_Software, 'used rounding approach - Mmtrx vs Stieltjes',RM_Type, 'which precisions we used',RM_Precision_list, 'what scaling set-up we used',RM_ScalingOptions };
% Info_ErrorCheckStuff = {'check cond. on norm',CheckErrMtrcs_norms, 'check cond. on entries of inv',CheckErrMtrcs_entrs, 'check cond. on P-regular splitting (eigvals)',CheckErrMtrcs_PregSplit};
% Info_PrblmAndDiscr = {'problem indices', ProblemChoices_list, '# interior grid points on the side of the square',list_of_nmb_int_gridcols};
% Output_Plotting = {'exact solution (only for the largest problem)', u_ExactSol, 'Error plots (only for the largest problem)', PlotData_ErrPlot, ...
%     'we saved all problem indicies and only the coarsest (=first) mesh resolution','indices of accuracies/digits for which we plot error', indsWhichPrecs_PlotErr, 'SM iteration indices for which we plot error', indsIter_PlotErr,};
% Output_Convergence = {'ConvCurvs for Schwarz methods', ConvCurves, 'approx ConvFacts for Schwarz methods', ConvFactApprox};
% Output_TheoryConds = {'condition ||F_imathcalAiInv|| < 1',ErrMtrcsNrms_CondStsfd, 'condition mathcalAiInv >= mathcalAiInvF_imathcalAiInv',ErrMtrcsEntrs_CondStsfd, 'condition lambdamin(A) >= lambdamax(Ei)',ErrMtrcsPregSplits_CondStsfd};
% 
% MyData = {'specs used for Schwarz method',Info_SchwarzMethodStuff, ... 
%     'specs used for rounding/multiprec',Info_RoundingStuff, ... 
%     'specs used for theory checks',Info_ErrorCheckStuff, ...
%     'specs for discretization and problem generation', Info_PrblmAndDiscr, ...
%     'outputs for plotting solutions and errors', Output_Plotting, ...
%     'outputs for convergence curves and convergence factors', Output_Convergence, ...
%     'outputs for checking theoretical conditions being satisfied', Output_TheoryConds};


clear; clc; close('all')


RM_Software = "chop";
%RM_Software = "advanpix";
Nmax = 55; RM_ScalingToRange = true; SM_RandomRHS = true;

% s_LoadString = "SavedData_mpSM_NonSymmPrblms_N" + Nmax + "_" +  RM_Software + ".mat"; LoadedData = load(s_LoadString).MyData; NonSymmPrblms = true;
% s_LoadString = "SavedData_mpSM_SymmPrblms_N" + Nmax + "_" +  RM_Software + ".mat"; LoadedData = load(s_LoadString).MyData; NonSymmPrblms = false;

% s_LoadString = "SavedData_mpSM_NonSymmPrblms_N" + Nmax + "_" +  RM_Software + "_Scaling" + RM_ScalingToRange + ".mat"; LoadedData = load(s_LoadString).MyData; NonSymmPrblms = true;
% s_LoadString = "SavedData_mpSM_SymmPrblms_N" + Nmax + "_" +  RM_Software + "_Scaling" + RM_ScalingToRange + ".mat"; LoadedData = load(s_LoadString).MyData; NonSymmPrblms = true;

s_LoadString = "SavedData_mpSM_NonSymmPrblms_N" + Nmax + "_" +  RM_Software + "_Scaling" + RM_ScalingToRange + "_RandomRHS" + SM_RandomRHS + ".mat"; LoadedData = load(s_LoadString).MyData; NonSymmPrblms = true;
s_LoadString = "SavedData_mpSM_SymmPrblms_N" + Nmax + "_" +  RM_Software + "_Scaling" + RM_ScalingToRange + "_RandomRHS" + SM_RandomRHS + ".mat"; LoadedData = load(s_LoadString).MyData; NonSymmPrblms = false;

 

Info_SchwarzMethodStuff = LoadedData{2};
SM_whichMethods = Info_SchwarzMethodStuff{2};

Info_RoundingStuff = LoadedData{4};
RM_whichPrecisions = Info_RoundingStuff{6};

Info_ErrorCheckStuff = LoadedData{6};
CheckErrMtrcs_norms = Info_ErrorCheckStuff{2}; CheckErrMtrcs_entrs = Info_ErrorCheckStuff{4}; CheckErrMtrcs_PregSplit = Info_ErrorCheckStuff{6};

Info_PrblmAndDiscr = LoadedData{8};
ProblemChoices_list = Info_PrblmAndDiscr{2}; list_of_nmb_int_gridcols = Info_PrblmAndDiscr{4};

Output_Plotting = LoadedData{10};
PlotData_ExactSol = Output_Plotting{2}; PlotData_ErrPlot = Output_Plotting{4}; indsWhichPrecs_PlotErr = Output_Plotting{7}; indsIter_PlotErr = Output_Plotting{9};

Output_Convergence = LoadedData{12};
ConvCurves = Output_Convergence{2}; ConvFactApprox = Output_Convergence{4};

Output_TheoryConds = LoadedData{14};
ErrMtrcsNrms_CondStsfd = Output_TheoryConds{2}; ErrMtrcsEntrs_CondStsfd = Output_TheoryConds{4}; ErrMtrcsPregSplits_CondStsfd = Output_TheoryConds{6};








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Plotting     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MyColors = {"#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F","#000000","#FF0000","#00FF00","#0000FF","#FFFF00","#FF00FF","#00FFFF"};
MyMarkers = {'d','^','o'}; MyLines = {'--','-.',':','-'};



%%% Plot the ConvFact as a function of number of digits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; t(1) = tiledlayout(2,3,'TileSpacing','compact','Padding','Compact');
ind_FrstRowWhichMethd = 3;
ind_MeshSizeToPlot = 1;
nmb_iters_ToPlot = 60; iters_mesh = 0:nmb_iters_ToPlot;
nmb_DigsToPlot = 3;

fontsize_ticks = 18; fontsize_axislabel = 25; fontsize_titles = 30;

%%% top row: convergence curves of MS for different d_{\ell} for the nonsymmetrical problems
for ind_whichPrblm = 1:length(ProblemChoices_list)
    curr_ax = nexttile(ind_whichPrblm); 
    for ind_nmbdig = nmb_DigsToPlot:-1:1
        PlotData = reshape(ConvCurves(ind_whichPrblm,ind_MeshSizeToPlot,ind_FrstRowWhichMethd,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
        semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker',MyMarkers{ind_FrstRowWhichMethd},'MarkerSize',12,'LineWidth',2); hold on;
    end
    if NonSymmPrblms, MyTitle = "Problem " + ind_whichPrblm; else, NewPrblmInd=3+ind_whichPrblm; MyTitle = "Problem " + NewPrblmInd; end
    title(MyTitle,'FontSize',fontsize_titles, 'interpreter', 'latex','FontWeight','bold'); 
    xlabel('iteration','FontSize',fontsize_ticks,'interpreter', 'latex'); %ylabel('2-norm of the error', 'FontSize',25,'interpreter', 'latex');
    set(curr_ax.XAxis,'FontSize',fontsize_ticks); set(curr_ax.XLabel,'FontSize',fontsize_axislabel);
    set(curr_ax.YAxis,'FontSize',fontsize_ticks); set(curr_ax.YLabel,'FontSize',fontsize_axislabel);
end


%%% bottom row: observed convergence factors of MS,dAS,RAS for different d_{\ell} for Problem 3,4,5
Precisions_mesh = 1:length(RM_whichPrecisions);



for ind_whichPrblm = 1:length(ProblemChoices_list)

    curr_ax = nexttile(3+ind_whichPrblm);

    for ind_whichMethd = 1:length(SM_whichMethods)
        PlotData = reshape(ConvFactApprox(ind_whichPrblm,ind_MeshSizeToPlot,ind_whichMethd,:) ,length(RM_whichPrecisions),1);
        plot( Precisions_mesh, PlotData, 'Color',[MyColors{7}],'Marker',MyMarkers{ind_whichMethd},'MarkerSize',14,'LineWidth',2); hold on;
        for ind_WhichPrec = 1:length(RM_whichPrecisions)

            if CheckErrMtrcs_norms && ~CheckErrMtrcs_PregSplit && ~CheckErrMtrcs_entrs
                if ErrMtrcsNrms_CondStsfd(ind_whichPrblm,ind_MeshSizeToPlot,ind_WhichPrec)
                    plot(ind_WhichPrec,ConvFactApprox(ind_whichPrblm,ind_MeshSizeToPlot,ind_whichMethd,ind_WhichPrec),'MarkerEdgeColor',[MyColors{7}], ...
                    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker',MyMarkers{ind_whichMethd},'MarkerSize',14,'LineWidth',2); hold on;
                end
            elseif CheckErrMtrcs_norms && CheckErrMtrcs_entrs && ~CheckErrMtrcs_PregSplit
                if ErrMtrcsNrms_CondStsfd(ind_whichPrblm,ind_MeshSizeToPlot,ind_WhichPrec) && ErrMtrcsEntrs_CondStsfd(ind_whichPrblm,ind_MeshSizeToPlot,ind_WhichPrec)
                    plot(ind_WhichPrec,ConvFactApprox(ind_whichPrblm,ind_MeshSizeToPlot,ind_whichMethd,ind_WhichPrec),'MarkerEdgeColor',[MyColors{7}], ...
                    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker',MyMarkers{ind_whichMethd},'MarkerSize',14,'LineWidth',2); hold on;
                end
            elseif CheckErrMtrcs_norms && ~CheckErrMtrcs_entrs && CheckErrMtrcs_PregSplit
                if ErrMtrcsNrms_CondStsfd(ind_whichPrblm,ind_MeshSizeToPlot,ind_WhichPrec) && ErrMtrcsPregSplits_CondStsfd(ind_whichPrblm,ind_MeshSizeToPlot,ind_WhichPrec)
                    plot(ind_WhichPrec,ConvFactApprox(ind_whichPrblm,ind_MeshSizeToPlot,ind_whichMethd,ind_WhichPrec),'MarkerEdgeColor',[MyColors{7}], ...
                    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker',MyMarkers{ind_whichMethd},'MarkerSize',14,'LineWidth',2); hold on;
                end
            elseif CheckErrMtrcs_norms && CheckErrMtrcs_entrs && CheckErrMtrcs_PregSplit
                if ErrMtrcsNrms_CondStsfd(ind_whichPrblm,ind_MeshSizeToPlot,ind_WhichPrec) && ErrMtrcsEntrs_CondStsfd(ind_whichPrblm,ind_MeshSizeToPlot,ind_WhichPrec) && ErrMtrcsPregSplits_CondStsfd(ind_whichPrblm,ind_MeshSizeToPlot,ind_WhichPrec)
                    plot(ind_WhichPrec,ConvFactApprox(ind_whichPrblm,ind_MeshSizeToPlot,ind_whichMethd,ind_WhichPrec),'MarkerEdgeColor',[MyColors{7}], ...
                    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker',MyMarkers{ind_whichMethd},'MarkerSize',14,'LineWidth',2); hold on;
                end

            end
        end

    end

    if strcmp(RM_Software,'chop')
        xticks(Precisions_mesh); xticklabels(RM_whichPrecisions);
    elseif strcmp(RM_Software,'advanpix')
        xlabel('$d_{\ell}$','FontSize',fontsize_axislabel, 'interpreter', 'latex');
    end    
    set(curr_ax.XAxis,'FontSize',fontsize_ticks); set(curr_ax.XLabel,'FontSize',fontsize_axislabel);
    set(curr_ax.YAxis,'FontSize',fontsize_ticks); set(curr_ax.YLabel,'FontSize',fontsize_axislabel);
end


%%% legend
if strcmp(RM_Software,'chop'), LegendLabels = cell(1,3); elseif strcmp(RM_Software,'advanpix'), LegendLabels = cell(1,3+nmb_DigsToPlot); end
LabelsForMethods = {' $\; \mathrm{dAS \; with} \; \theta = \frac{1}{3} \quad$',' $\; \mathrm{RAS} \quad$','$\; \mathrm{MS} \quad$'};

for ind_whichMethd = 1:length(SM_whichMethods)
plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker',MyMarkers{ind_whichMethd},'MarkerSize',16,'LineWidth',2); hold on;
LegendLabels{ind_whichMethd} = LabelsForMethods{ind_whichMethd}; LegendHandles(ind_whichMethd) = plt1;
end

for ind_WhichPrec = 1:min([length(RM_whichPrecisions),nmb_DigsToPlot])
    plt1 = plot(nan,nan,'Color',[MyColors{ind_WhichPrec}],'LineWidth',6); hold on;
    LegendHandles(3+ind_WhichPrec) = plt1; 
    if strcmp(RM_Software,'chop')
        LegendLabels{3+ind_WhichPrec} = RM_whichPrecisions{ind_WhichPrec};
    elseif strcmp(RM_Software,'advanpix')
        LegendLabels{3+ind_WhichPrec} = "$d_{\ell}$ = " + RM_whichPrecisions(ind_WhichPrec);
    end
end

set(groot,'defaultLegendInterpreter','latex');
Lgnd = legend(LegendHandles,LegendLabels,'Location','northoutside','Orientation','Horizontal'); Lgnd.NumColumns = 10; Lgnd.Layout.Tile = 'North';
fontsize(Lgnd,fontsize_axislabel,'points');






%%% Plot the errors and their evolution over iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_whichProblmToPlot = 2;
ind_whichMthdToPlot = 3;
ind_MeshSizeToPlot = 1;
indsIter_PlotErr = indsIter_PlotErr(1:3);
nmb_int_gridcols = list_of_nmb_int_gridcols(ind_MeshSizeToPlot); h = 1/(nmb_int_gridcols+1); u_plot = zeros(nmb_int_gridcols+2); x_mesh = 0:h:1; angle1 = 224; angle2 = 35;

nmbRowsTiles = length(indsWhichPrecs_PlotErr); nmbColsTiles = length(indsIter_PlotErr);
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');
fontsize_ticks = 18; fontsize_axislabel = 25; fontsize_titles = 30;


for ind_prec = 1:length(indsWhichPrecs_PlotErr)
    for ind_iter = 1:length(indsIter_PlotErr)

        curr_ax = nexttile();
        u_err = reshape( PlotData_ErrPlot{ind_whichProblmToPlot,ind_whichMthdToPlot,ind_prec,ind_iter}, 1,nmb_int_gridcols^2);
        PlotData = u_plot;
        PlotData(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = reshape(u_err,nmb_int_gridcols,nmb_int_gridcols)'; 
        mesh(x_mesh,x_mesh,PlotData); view(angle1,angle2); 
        max_err = max(max(abs(u_err))); max_err_exp = log10(max_err); if ind_prec == 2, max_err_exp = max_err_exp - 1; end
        ax = gca; ax.ZAxis.Exponent = floor(max_err_exp);

        if ind_prec == 1 
            mylabel = "iteration " + indsIter_PlotErr(ind_iter); title(mylabel, 'interpreter', 'latex','FontSize',fontsize_titles);
        end
        if ind_prec == length(indsWhichPrecs_PlotErr)
            xlabel('$x_1$','interpreter', 'latex'); ylabel('$x_2$','interpreter', 'latex');
        else
            xticklabels({}); yticklabels({});
        end
        set(curr_ax.XAxis,'FontSize',fontsize_ticks); set(curr_ax.XLabel,'FontSize',fontsize_ticks);
        set(curr_ax.YAxis,'FontSize',fontsize_ticks); set(curr_ax.YLabel,'FontSize',fontsize_ticks);
        set(curr_ax.ZAxis,'FontSize',fontsize_ticks); set(curr_ax.ZLabel,'FontSize',fontsize_ticks);

        % if ind_iter == 1, mylabel = append('$d_{\ell}$ = ',num2str(indsDigs_PlotErr_Prblm3(ind_dig))); zlabel(mylabel, 'interpreter', 'latex','FontSize',26); 
        %    zl=get(gca,'zlabel'); pzl = get(zl,'position'); pzl(1) = 1.001*pzl(1); set(zl,'position',pzl);
        % end

        u_plot = zeros(nmb_int_gridcols+2);

    end
end



%%% Plot the exact solution (largest grid, last problem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_whichProblmToPlot = 3;
nmb_int_gridcols = list_of_nmb_int_gridcols(1); h = 1/(nmb_int_gridcols+1); u_plot = zeros(nmb_int_gridcols+2); x_mesh = 0:h:1; angle1 = 224; angle2 = 35;

figure; fontsize_ticks = 18; fontsize_axislabel = 25; fontsize_titles = 30;

curr_ax = gca;
u_sol = reshape( PlotData_ExactSol{ind_whichProblmToPlot}, 1,nmb_int_gridcols^2);
PlotData = u_plot;
PlotData(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = reshape(u_sol,nmb_int_gridcols,nmb_int_gridcols)'; 
mesh(x_mesh,x_mesh,PlotData); view(angle1,angle2); 
max_err = max(max(abs(u_sol))); max_err_exp = log10(max_err); if ind_prec == 2, max_err_exp = max_err_exp - 1; end
ax = gca; ax.ZAxis.Exponent = floor(max_err_exp);

xlabel('$x_1$','interpreter', 'latex'); ylabel('$x_2$','interpreter', 'latex');
set(curr_ax.XAxis,'FontSize',fontsize_ticks); set(curr_ax.XLabel,'FontSize',fontsize_ticks);
set(curr_ax.YAxis,'FontSize',fontsize_ticks); set(curr_ax.YLabel,'FontSize',fontsize_ticks);
set(curr_ax.ZAxis,'FontSize',fontsize_ticks); set(curr_ax.ZLabel,'FontSize',fontsize_ticks);














%%% plot interaction of "d_{\ell} " and "h"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontsize_ticks = 18; fontsize_axislabel = 25; fontsize_titles = 30;

nmbRowsTiles = 3; nmbColsTiles = 3; 
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');
LegendHandles = []; LegendLabels = {};

angle1 = 227; angle2 = 26;
% Cbar_Ticks = []; Cbar_TickLabels = cell(length(Cbar_Ticks));
Cbar_Ticks = [1,3,5,7,9,11,13,15]; 
Cbar_TickLabels = cell(length(Cbar_Ticks),1);
for ind = 1:length(Cbar_Ticks)
    Cbar_TickLabels{ind} = int2str( list_of_nmb_int_gridcols(Cbar_Ticks(ind))^2 );
end


% Observed convergence factor for diff d_{\ell} and nmb_int_gridcols for Problems 1,2,3 & dAS,RAS,MS 
if CheckErrMtrcs_norms && ~CheckErrMtrcs_entrs && ~CheckErrMtrcs_PregSplit
    [ErrMtrcsNrms_MyFrstIndcs] = FrstIndxToSatisfyCond(ErrMtrcsNrms_CondStsfd,nan,nan); %%% remember "ErrMtrcsNrms_CondStsfd(ind_whichPrblm,ind_MeshSizeToPlot,ind_WhichPrec)"
elseif CheckErrMtrcs_norms && ~CheckErrMtrcs_entrs && CheckErrMtrcs_PregSplit
    [ErrMtrcsNrms_MyFrstIndcs] = FrstIndxToSatisfyCond(ErrMtrcsNrms_CondStsfd,nan,nan); %%% remember "ErrMtrcsNrms_CondStsfd(ind_whichPrblm,ind_MeshSizeToPlot,ind_WhichPrec)"
    [ErrMtrcsPregSplit_MyFrstIndcs] = FrstIndxToSatisfyCond(ErrMtrcsPregSplits_CondStsfd,nan,nan); %%% remember "ErrMtrcsNrms_CondStsfd(ind_whichPrblm,ind_MeshSizeToPlot,ind_WhichPrec)"
end

for ind_whichPrblm = 1:length(ProblemChoices_list)
    for ind_whichMethd = 1:length(SM_whichMethods)

        curr_ax = nexttile( (ind_whichPrblm-1)*3 + ind_whichMethd );

        PlotData = reshape(ConvFactApprox(ind_whichPrblm,:,ind_whichMethd,:) ,length(list_of_nmb_int_gridcols),length(RM_whichPrecisions));
        ribbon(PlotData'); hold on; 

        for ind_MeshSizeToPlot = 1:length(list_of_nmb_int_gridcols)

            if CheckErrMtrcs_norms && ~CheckErrMtrcs_entrs && ~CheckErrMtrcs_PregSplit
                x = ind_MeshSizeToPlot; y = ErrMtrcsNrms_MyFrstIndcs(ind_whichPrblm,ind_MeshSizeToPlot); z = ConvFactApprox(ind_whichPrblm,x,ind_whichMethd,y);
                plt = plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.7 .7 .7]'); hold on;
                if ind_MeshSizeToPlot == 1 && ind_whichPrblm == 1 && ind_whichMethd == 1, LegendHandles(1) = plt; LegendLabels{1} = 'first $d_{\ell}$ such that $\| \mathcal{A}_i^{-1} \mathtt{F}_i \| < 1$'; end
            elseif CheckErrMtrcs_norms && ~CheckErrMtrcs_entrs && CheckErrMtrcs_PregSplit
                x = ind_MeshSizeToPlot; y = ErrMtrcsNrms_MyFrstIndcs(ind_whichPrblm,ind_MeshSizeToPlot); z = ConvFactApprox(ind_whichPrblm,x,ind_whichMethd,y);
                plt = plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.7 .7 .7]'); hold on;
                if ind_MeshSizeToPlot == 1 && ind_whichPrblm == 1 && ind_whichMethd == 1, LegendHandles(1) = plt; LegendLabels{1} = 'first $d_{\ell}$ such that $\| \mathcal{A}_i^{-1} \mathtt{F}_i \| < 1$'; end
                x = ind_MeshSizeToPlot; y = ErrMtrcsPregSplit_MyFrstIndcs(ind_whichPrblm,ind_MeshSizeToPlot); z = ConvFactApprox(ind_whichPrblm,x,ind_whichMethd,y);
                plt = plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.95 .95 .95]'); hold on;
                if ind_MeshSizeToPlot == 1 && ind_whichPrblm == 1 && ind_whichMethd == 1, LegendHandles(2) = plt; LegendLabels{2} = 'first $d_{\ell}$ such that $\lambda_{min}( \mathcal{A}_i ) \geq 2 \left| \lambda_{-\infty}( \mathtt{F}_i ) \right|$'; end
            end
        end
        view(angle1,angle2); 
        if ind_whichPrblm == 1 && ind_whichMethd == 1, title('dAS','FontSize',fontsize_titles, 'interpreter', 'latex');
        elseif ind_whichPrblm == 1 && ind_whichMethd == 2, title('RAS','FontSize',fontsize_titles, 'interpreter', 'latex');
        elseif ind_whichPrblm == 1 && ind_whichMethd == 3, title('MS','FontSize',fontsize_titles, 'interpreter', 'latex'); end
        xticks([]); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({});

        if ind_whichPrblm == length(ProblemChoices_list)
            yticklabels({'1','4','8','12','16'}); ylabel('$d_{\ell}$','FontSize',24,'interpreter', 'latex'); ytickangle(curr_ax,0)
        end
        
        if NonSymmPrblms
            if ind_whichPrblm == 1 && ind_whichMethd == 1, zlabel('Problem 1','FontSize',24,'interpreter', 'latex');
            elseif ind_whichPrblm == 2 && ind_whichMethd == 1, zlabel('Problem 2','FontSize',24,'interpreter', 'latex');
            elseif ind_whichPrblm == 3 && ind_whichMethd == 1, zlabel('Problem 3','FontSize',24,'interpreter', 'latex'); end
        else
            if ind_whichPrblm == 1 && ind_whichMethd == 1, zlabel('Problem 4','FontSize',24,'interpreter', 'latex');
            elseif ind_whichPrblm == 2 && ind_whichMethd == 1, zlabel('Problem 5','FontSize',24,'interpreter', 'latex');
            elseif ind_whichPrblm == 3 && ind_whichMethd == 1, zlabel('Problem 6','FontSize',24,'interpreter', 'latex'); end
        end
        set(curr_ax.YAxis,'FontSize',fontsize_ticks); set(curr_ax.YLabel,'FontSize',fontsize_axislabel);
        set(curr_ax.ZAxis,'FontSize',fontsize_ticks); set(curr_ax.ZLabel,'FontSize',fontsize_axislabel);

    end
end



% colorbar
cbar = colorbar('Location','eastoutside'); cbar.Label.String= '$N$'; cbar.Label.Color= 'black'; 
cbar.Ticks = Cbar_Ticks; cbar.TickLabels = Cbar_TickLabels; 
cbar.Label.Interpreter = 'latex'; cbar.TickLabelInterpreter = 'latex';
cbar.Layout.Tile = 'east';
set(cbar,'FontSize',fontsize_ticks);

% Legend
set(groot,'defaultLegendInterpreter','latex');
Lgnd = legend(LegendHandles,LegendLabels,'Location','northoutside','Orientation','Horizontal'); Lgnd.NumColumns = 7; Lgnd.Layout.Tile = 'North';
fontsize(Lgnd,fontsize_axislabel,'points'); 
















function[MyFrstIndcs] = FrstIndxToSatisfyCond(ArrayOfConds_1,ArrayOfConds_2,ArrayOfConds_3)
[d1,d2,~] = size(ArrayOfConds_1); MyFrstIndcs = nan(d1,d2);
if any(isnan(ArrayOfConds_3))

    if any(isnan(ArrayOfConds_2))

        if any(isnan(ArrayOfConds_1))
            MyFrstIndcs = nan;
        else
            for ind1=1:d1
                for ind2=1:d2
                    MyFrstInd = find(ArrayOfConds_1(ind1,ind2,:),1,'first');
                    if isempty(MyFrstInd), MyFrstIndcs(ind1,ind2) = nan; else, MyFrstIndcs(ind1,ind2) = MyFrstInd; end
                end
            end
        end

    else
        for ind1=1:d1
            for ind2=1:d2
                
                MyFrstInd_Conds1 = find(ArrayOfConds_1(ind1,ind2,:),1,'first'); MyFrstInd_Conds2 = find(ArrayOfConds_2(ind1,ind2,:),1,'first');
                if isempty(MyFrstInd_Conds1) || isempty(MyFrstInd_Conds2)
                    MyFrstIndcs(ind1,ind2) = nan;
                else
                    MyFrstIndcs(ind1,ind2) = max([MyFrstInd_Conds1,MyFrstInd_Conds2]);
                end
                
            end
        end

    end

else
    for ind1=1:d1
        for ind2=1:d2
            
                MyFrstInd_Conds1 = find(ArrayOfConds_1(ind1,ind2,:),1,'first'); MyFrstInd_Conds2 = find(ArrayOfConds_2(ind1,ind2,:),1,'first'); MyFrstInd_Conds3 = find(ArrayOfConds_3(ind1,ind2,:),1,'first');
                if isempty(MyFrstInd_Conds1) || isempty(MyFrstInd_Conds2) || isempty(MyFrstInd_Conds3)
                    MyFrstIndcs(ind1,ind2) = nan;
                else
                    MyFrstIndcs(ind1,ind2) = max([MyFrstInd_Conds1,MyFrstInd_Conds2,MyFrstInd_Conds3]);
                end

        end
    end

end
end