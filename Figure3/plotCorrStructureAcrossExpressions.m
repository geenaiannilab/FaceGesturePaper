
%%%%%%%%
%%%%%%%% PLOTTING Figure 4C, D
%%%%%%%%  Kinematic decoding and non-preserved neural correlations
%%%%%%%%
%%%%%%%% Plots neural signal correlation between all neuron pairs over time, per behavior (trial-avgeraged) 
%%%%%%%%    Each matrix of cell-cell correlations is clustered, either by 
%%%%%%%%     hierarchical clustering of its own structure (top row), or a
%%%%%%%%     unified clustering defined by the lipsmack matrix (bottom row)
%%%%%%%%     
%%%%%%%  Plots all R^2 values (1 per day, between pairwise correlations:
%%%%%%        threat-lipsmack, threat-chew, lipsmack-chew); 
%%%%%%%%% GRI 09/11/2025 

close all; clear all 

set(0,'defaultAxesFontSize', 24); % bc im blind 
set(0,'defaultAxesFontWeight', 'bold'); 

date = '171128';
subject ='Thor';
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/CorrMatricesV2']);
regions = {'All'};
bhvs2plot = [1 2 4];
bhvStrings = {'Threat','Lipsmack','Chew'};

% region colormap
cmap1 = [0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250
    0.5 0.5 0.5];

% for each region, plot the self-sorted & unified sorting signal
% correlations, per behavior 
%%
for rr = 1:length(regions)
    data = load([workdir '/' regions{rr} '.mat']);
    
    taxis = data.taxis;

     % plot pairwise correlations
    r2_1v2 = rsquare(data.prCorr{:,1}', data.prCorr{:,2}');
    r2_1v3 = rsquare(data.prCorr{:,1}', data.prCorr{:,3}');
    r2_2v3 = rsquare(data.prCorr{:,2}', data.prCorr{:,3}');

    disp(['Region is ' regions{rr}])
    disp([regions{rr} ': Mean R2 = ' num2str(mean([r2_1v2 r2_1v3 r2_2v3]))])

    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    fig.Position = [0    0.4399    1.0000    0.4409];

    subplot 131; scatter(data.prCorr{:,1},data.prCorr{:,2},'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
    xlabel('Pairwise Correlation (Threat)','FontSize',20); ylabel('Pairwise Correlation (Lipsmack)','FontSize',20); xlim([-1 1]); ylim([-1 1 ]);
    axes = gca; axes.XTick = [-1 0 1]; axes.YTick = [-1 0 1];
    %text(-1,-1,['R2 = ' num2str( r2_1v2 )]);
    axes.FontWeight = 'bold';
    disp(['R2 = ' num2str( r2_1v2 )]);

    subplot 132; scatter(data.prCorr{:,1},data.prCorr{:,3},'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
    xlabel('Pairwise Correlation (Threat)','FontSize',20); ylabel('Pairwise Correlation (Chew)','FontSize',20);xlim([-1 1]); ylim([-1 1 ]);
    axes = gca; axes.XTick = [-1 0 1]; axes.YTick = [-1 0 1];
    %text(-1,-1,['R2 = ' num2str( r2_1v3 )]);
    axes.FontWeight = 'bold';
    disp(['R2 = ' num2str( r2_1v3 )]);

    subplot 133; scatter(data.prCorr{:,2},data.prCorr{:,3},'MarkerFaceColor','k','MarkerEdgeColor','k',...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    xlabel('Pairwise Correlation (Lipsmack)','FontSize',20); ylabel('Pairwise Correlation (Chew)','FontSize',20);xlim([-1 1]); ylim([-1 1 ]);
    axes = gca; axes.XTick = [-1 0 1]; axes.YTick = [-1 0 1];
    axes.FontWeight = 'bold';
    %text(-1,-1,['R2 = ' num2str( r2_2v3 )]);
    disp(['R2 = ' num2str( r2_2v3 )]);

    sgtitle(['Neuron Pairwise Correlations: ' regions{rr}],'FontSize',36,'FontWeight','bold');

    %plotFig = gcf;
    %savefig(plotFig, ['~/Dropbox/PhD/Manuscript/subFigurePDFs/' regions{rr} '_' date '_neuralSignalCorr_scatter.fig']);
    %exportgraphics(plotFig, ['~/Dropbox/PhD/Manuscript/subFigurePDFs/' regions{rr} '_' date '_neuralSignalCorr_scatter.pdf'], 'ContentType', 'vector');

end

% correlation matrices 
for rr = 1:length(regions)
    
    data = load([workdir '/' regions{rr} '.mat']);

   clear axes;

    fig = figure; fig.Position =  [ 476          97        1037         769];
    ha = tight_subplot(2,3,[.05 .05], [0.05 0.1]);

    for bhv = 1:length(bhvs2plot)

        axes(ha(bhv));
        
        % self-defined clusters 
      
        imagesc(data.signalCorrOverTime(bhv).selfSort.data); 
        ax1 = gca;
        ax1.XTick = []; %1:length(spikeLabels);
        ax1.YTick = []; %1:length(spikeLabels);

        colorbarpzn(0,2,'full',1,'off')
        xlabel('Cell ID','FontSize',20); ylabel('Cell ID','FontSize',20)
        title([bhvStrings{bhv}],'FontSize',28)

        % uniformly-sorted clusters 
        axes(ha(bhv+3));
        imagesc(data.signalCorrOverTime(bhv).unifiedSort.data); 
        ax1 = gca;
        ax1.XTick = [];%;1:length(spikeLabels);
        ax1.YTick = [];%1:length(spikeLabels);
        colorbarpzn(0,2,'full',1,'off')
        xlabel('Cell ID','FontSize',20); ylabel('Cell ID','FontSize',20)

        t = sgtitle(regions{rr},'FontSize',36,'FontWeight','bold');

    end % end bhvs 

    c = colorbar;
    c.Position = [0.96 0.09 0.01 0.8];
    c.Ticks = [0 1 2];
    c.FontWeight= 'bold';

    %savefig(fig, ['~/Dropbox/PhD/Manuscript/subFigurePDFs/' regions{rr} '_' date '_neuralSignalCorr.fig']);
    %exportgraphics(fig, ['~/Dropbox/PhD/Manuscript/subFigurePDFs/' regions{rr} '_' date '_neuralSignalCorr.pdf'], 'ContentType', 'vector');
end % end regions
%%

% Plot the R^2 
thorData = load('/Users/geena/Dropbox/PhD/SUAinfo/Thor_R2ValsAcrossDays.mat');
barneyData = load('/Users/geena/Dropbox/PhD/SUAinfo/Barney_R2ValsAcrossDays.mat');

thorMeans_1v2 = mean(thorData.r2_1v2,1); 
thorMeans_1v3 = mean(thorData.r2_1v3,1);
thorMeans_2v3 = mean(thorData.r2_2v3,1);

barneyMeans_1v2 = mean(barneyData.r2_1v2,1);
barneyMeans_1v3 = mean(barneyData.r2_1v3,1);
barneyMeans_2v3 = mean(barneyData.r2_2v3,1);

combinedData_1v2 = [thorData.r2_1v2; barneyData.r2_1v2];
combinedData_1v3 = [thorData.r2_1v3; barneyData.r2_1v3];
combinedData_2v3 = [thorData.r2_2v3; barneyData.r2_2v3];

combinedPval_1v2 = [thorData.p_1v2; barneyData.p_1v2];
combinedPval_1v3 = [thorData.p_1v3; barneyData.p_1v3];
combinedPval_2v3 = [thorData.p_2v3; barneyData.p_2v3];

combinedMeans_1v2 = mean([thorData.r2_1v2 ; barneyData.r2_1v2],1);
combinedMeans_1v3 = mean([thorData.r2_1v3 ; barneyData.r2_1v3],1);
combinedMeans_2v3 = mean([thorData.r2_2v3 ; barneyData.r2_2v3],1);

combinedStd_1v2 = std([thorData.r2_1v2 ; barneyData.r2_1v2],[],1);
combinedStd_1v3 = std([thorData.r2_1v3 ; barneyData.r2_1v3],[],1);
combinedStd_2v3 = std([thorData.r2_2v3 ; barneyData.r2_2v3],[],1);

combinedMeansAllBhv = mean([thorData.r2_1v2 ; thorData.r2_1v3; thorData.r2_2v3; barneyData.r2_1v2 ;barneyData.r2_1v3; barneyData.r2_2v3 ],1);
combinedStdAllBhv = std([thorData.r2_1v2 ; thorData.r2_1v3; thorData.r2_2v3; barneyData.r2_1v2 ;barneyData.r2_1v3; barneyData.r2_2v3 ],[],1);

combinedMeansAllBhvAllRegions = mean([thorData.r2_1v2 ; thorData.r2_1v3; thorData.r2_2v3; barneyData.r2_1v2 ;barneyData.r2_1v3; barneyData.r2_2v3 ],'all');
combinedStdAllBhvAllRegions = std([thorData.r2_1v2 ; thorData.r2_1v3; thorData.r2_2v3; barneyData.r2_1v2 ;barneyData.r2_1v3; barneyData.r2_2v3 ],[],'all');


% Plot R2 for threat v lipsmack (bhv1v2)
figure; 
counter = 1;

for dd = 1:size(combinedData_1v2,1)
    for rr = 1:length(regions)
        
        
        jitter = rand(1,1)*0.5;
        if rem(counter, 2) == 0
            x = rr*2 + jitter;
        else
             x = rr*2 - jitter;
        end
        s = swarmchart(x, combinedData_1v2(dd,rr),250,cmap1(rr,:),'filled'); hold on;
        alpha = (1-combinedPval_1v2(dd,rr))+0.1;
        if alpha > 1
            s.MarkerFaceAlpha = 1;
        else
            s.MarkerFaceAlpha = alpha;
        end
        
        counter = counter+1;

    end
end

% Plot R2 for threat V chew (bhv1v3)
counter = 1;

for dd = 1:size(combinedData_1v3,1)
    for rr = 1:length(regions)
        
        
        jitter = rand(1,1)*0.5;
        if rem(counter, 2) == 0
            x = rr*2 + jitter;
        else
             x = rr*2 - jitter;
        end
        s = swarmchart(x, combinedData_1v3(dd,rr),250,cmap1(rr,:),'filled'); hold on;
        alpha = (1-combinedPval_1v3(dd,rr))+0.1;
        if alpha > 1
            s.MarkerFaceAlpha = 1;
        else
            s.MarkerFaceAlpha = alpha;
        end
        
        counter = counter+1;

    end
end

% Plot R2 for lipsmack V chew (bhv2v3)

counter = 1;

for dd = 1:size(combinedData_2v3,1)
    for rr = 1:length(regions)
        
        
        jitter = rand(1,1)*0.5;
        if rem(counter, 2) == 0
            x = rr*2 + jitter;
        else
             x = rr*2 - jitter;
        end
        s = swarmchart(x, combinedData_2v3(dd,rr),250,cmap1(rr,:),'filled'); hold on;
        alpha = (1-combinedPval_2v3(dd,rr))+0.1;
        if alpha > 1
            s.MarkerFaceAlpha = 1;
        else
            s.MarkerFaceAlpha = alpha;
        end
        
        counter = counter+1;

    end
end

for rr = 1:length(regions)
    scatter(rr*2, combinedMeansAllBhv(rr),500,cmap1(rr,:),'+','linew',4);
    hold on 
end
hold off
ylabel('R-Squared');
title('Pairwise Neuron Correlations');

hold off; 
axes = gca;
axes.XTick = 2:2:10;
axes.XTickLabel = regions;

%%

