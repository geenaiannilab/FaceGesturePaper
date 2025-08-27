% written GI 
%%%     
%%% PUPRPOSE: plot the 1) signal correlation between all neuron pairs over time, per behavior (trial-avgeraged) 
%%%                  2) noise correlations between all neuron pairs, per behavior (per trial) 

close all; clear all 

set(0,'defaultAxesFontSize', 24); % bc im blind 
set(0,'defaultAxesFontWeight', 'bold'); 

date = '210704';
subject ='Barney';
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/CorrMatricesV2']);
regions = {'S1','M1','PMv','M3','All'};
bhvs2plot = [1 2 4];
bhvStrings = {'Threat','Lipsmack','Chew'};


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

    plotFig = gcf;
    savefig(plotFig, ['~/Dropbox/PhD/Manuscript/subFigurePDFs/' regions{rr} '_' date '_neuralSignalCorr_scatter.fig']);
    exportgraphics(plotFig, ['~/Dropbox/PhD/Manuscript/subFigurePDFs/' regions{rr} '_' date '_neuralSignalCorr_scatter.pdf'], 'ContentType', 'vector');

end


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

    savefig(fig, ['~/Dropbox/PhD/Manuscript/subFigurePDFs/' regions{rr} '_' date '_neuralSignalCorr.fig']);
    exportgraphics(fig, ['~/Dropbox/PhD/Manuscript/subFigurePDFs/' regions{rr} '_' date '_neuralSignalCorr.pdf'], 'ContentType', 'vector');
end % end regions
%%
