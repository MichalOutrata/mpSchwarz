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
RM_Software = "advanpix";
Nmax = 330; RM_ScalingToRange = true; SM_RandomRHS = true;

%s_LoadString = "SavedData_mpSMprec_NonSymmPrblms_N" + Nmax + "_" +  RM_Software + "_Scaling" + RM_ScalingToRange + "_RandomRHS" + SM_RandomRHS + ".mat"; LoadedData = load(s_LoadString).MyData; NonSymmPrblms = true;
s_LoadString = "SavedData_mpSMprec_SymmPrblms_N" + Nmax + "_" +  RM_Software + "_Scaling" + RM_ScalingToRange + "_RandomRHS" + SM_RandomRHS + ".mat"; LoadedData = load(s_LoadString).MyData; NonSymmPrblms = false;


Info_SchwarzMethodStuff = LoadedData{2};
SM_whichMethods = Info_SchwarzMethodStuff{2};

Info_GMRESStuff = LoadedData{4};
GMRES_nmbiter = Info_GMRESStuff{4};

Info_RoundingStuff = LoadedData{6};
RM_whichPrecisions = Info_RoundingStuff{6};

Info_PrblmAndDiscr = LoadedData{8};
ProblemChoices_list = Info_PrblmAndDiscr{2}; list_of_nmb_int_gridcols = Info_PrblmAndDiscr{4};

Output_Convergence = LoadedData{10};
GMRESprec_ConvCurves = Output_Convergence{2}; GMRESprec_NmbItToConv = Output_Convergence{4};









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Plotting     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MyColors = {"#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F","#000000","#FF0000","#00FF00","#0000FF","#FFFF00","#FF00FF","#00FFFF"};
MyMarkers = {'d','^','o'}; MyLines = {'--','-.',':','-'};



%%% Plot the ConvFact as a function of number of digits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; t(1) = tiledlayout(2,3,'TileSpacing','compact','Padding','Compact');
ind_whichMethod = 3;
ind_MeshSizeToPlot = 1;
nmb_iters_ToPlot = GMRES_nmbiter; 
iters_mesh = 0:nmb_iters_ToPlot;
nmb_DigsToPlot = 6; angle1 = 227; angle2 = 26;

Cbar_Ticks = [1,3,5,7,9,11,13,15]; Cbar_Ticks = 1:2:length(list_of_nmb_int_gridcols); 
Cbar_TickLabels = cell(length(Cbar_Ticks),1);
for ind = 1:length(Cbar_Ticks)
    Cbar_TickLabels{ind} = int2str( list_of_nmb_int_gridcols(Cbar_Ticks(ind))^2 );
end


fontsize_ticks = 18; fontsize_axislabel = 25; fontsize_titles = 30;

%%% top row: convergence curves of MS for different d_{\ell} for the nonsymmetrical problems
for ind_whichPrblm = 1:length(ProblemChoices_list)
    curr_ax = nexttile(ind_whichPrblm); 
    for ind_nmbdig = nmb_DigsToPlot:-1:1
        PlotData = reshape( GMRESprec_ConvCurves(ind_whichPrblm, ind_MeshSizeToPlot, ind_nmbdig, ind_whichMethod,1:nmb_iters_ToPlot+1), nmb_iters_ToPlot+1,1);
        semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker',MyMarkers{ind_whichMethod},'MarkerSize',12,'LineWidth',2); hold on;
    end
    if NonSymmPrblms, MyTitle = "Problem " + ind_whichPrblm; else, NewPrblmInd=3+ind_whichPrblm; MyTitle = "Problem " + NewPrblmInd; end
    title(MyTitle,'FontSize',fontsize_titles, 'interpreter', 'latex','FontWeight','bold'); 
    xlabel('iteration','FontSize',fontsize_ticks,'interpreter', 'latex'); %ylabel('2-norm of the error', 'FontSize',25,'interpreter', 'latex');
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


%%% bottom row: observed convergence factors of MS,dAS,RAS for different d_{\ell} for Problem 3,4,5
Precisions_mesh = 1:length(RM_whichPrecisions);


for ind_whichPrblm = 1:length(ProblemChoices_list)

    curr_ax = nexttile(3+ind_whichPrblm);
    PlotData = reshape(GMRESprec_NmbItToConv(ind_whichPrblm, :, :, ind_whichMethod) ,length(list_of_nmb_int_gridcols),length(RM_whichPrecisions));
    ribbon(PlotData'); view(angle1,angle2);


    if strcmp(RM_Software,'chop')
        xticks([]); yticks(Precisions_mesh); yticklabels(RM_whichPrecisions);
    elseif strcmp(RM_Software,'advanpix')
        xticks([]); ylabel('$d_{\ell}$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'});
    end
    if ind_whichPrblm==1, zlabel('\# GMRES iterations','FontSize',24,'interpreter', 'latex'); end
    set(curr_ax.XAxis,'FontSize',fontsize_ticks); set(curr_ax.XLabel,'FontSize',fontsize_axislabel);
    set(curr_ax.YAxis,'FontSize',fontsize_ticks); set(curr_ax.YLabel,'FontSize',fontsize_axislabel);
end
cbar = colorbar(); cbar.Label.String= '$N$'; cbar.Label.Color= 'black'; 
cbar.Ticks = Cbar_Ticks; cbar.TickLabels = Cbar_TickLabels; 
cbar.Label.Interpreter = 'latex'; cbar.TickLabelInterpreter = 'latex';
% cbar.Layout.Tile = 'east';


















%%% plot interaction of "d_{\ell} " and "h"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontsize_ticks = 18; fontsize_axislabel = 25; fontsize_titles = 30;

nmbRowsTiles = 3; nmbColsTiles = 3; 
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');
LegendHandles = []; LegendLabels = {};

for ind_whichPrblm = 1:length(ProblemChoices_list)
    for ind_whichMethd = 1:length(SM_whichMethods)

        curr_ax = nexttile( (ind_whichPrblm-1)*3 + ind_whichMethd );

        PlotData = reshape(GMRESprec_NmbItToConv(ind_whichPrblm, :, :, ind_whichMethd) ,length(list_of_nmb_int_gridcols),length(RM_whichPrecisions));
        ribbon(PlotData'); hold on; 

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

% % Legend
% set(groot,'defaultLegendInterpreter','latex');
% Lgnd = legend(LegendHandles,LegendLabels,'Location','northoutside','Orientation','Horizontal'); Lgnd.NumColumns = 7; Lgnd.Layout.Tile = 'North';
% fontsize(Lgnd,fontsize_axislabel,'points'); 











%%% plot interaction of "d_{\ell} " and GMRES convergence curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_Prblm = 2; % choose
ind_meshsize = length(list_of_nmb_int_gridcols); % choose
ind_whichMethod = 2; % choose

ClsstIntSqrt = floor(sqrt(length(RM_whichPrecisions))); nmbColsTiles = floor( length(RM_whichPrecisions) / ClsstIntSqrt );
if length(RM_whichPrecisions) <= 3, nmbColsTiles = length(RM_whichPrecisions); nmbRowsTiles = 1; else
    if floor( length(RM_whichPrecisions) / nmbColsTiles ) == length(RM_whichPrecisions) / nmbColsTiles, nmbRowsTiles = floor( length(RM_whichPrecisions) / nmbColsTiles ); 
    else, nmbRowsTiles = floor( length(RM_whichPrecisions) / nmbColsTiles ) + 1; end
end
figure(96); t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');
for ind_nmbdig = 1:length(RM_whichPrecisions)
    ax(ind_nmbdig) = nexttile();
    PlotData = reshape( GMRESprec_ConvCurves(ind_Prblm, ind_meshsize, ind_nmbdig, ind_whichMethod,:), GMRES_nmbiter+1,1);
    semilogy(PlotData,'ko-');
end
title(t, append('# grid points on one side = ',num2str(list_of_nmb_int_gridcols(ind_meshsize))));

% figure(97)
% plot(RM_Precision_list, reshape( GMRESprec_NmbItToConv(ind_Prblm, ind_meshsize, :, ind_whichMethod), size(RM_Precision_list) ),'ro-'); xlabel('# digits');ylabel('# GMRES it to conv'); title( append('# grid points on one side = ',num2str(list_of_nmb_int_gridcols(ind_meshsize))));