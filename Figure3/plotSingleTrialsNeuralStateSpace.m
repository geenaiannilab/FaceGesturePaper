
%%%%%%%% PLOTTING Figure 3A
%%%%%%%%  Kinematic decoding and non-preserved neural correlations
%%%%%%%%
%%%%%%%%  Plots all Spearman's rho values (large for preserved rank-order, small if not) 

% This script will calculate the 1) loadings of each neuron onto each PC 
%                                2) the projections of the data onto the
%                                PCs
%                                3) variance explained by each PC 

clear all ; 
data = load('/MatFiles/Fig3A.mat');
sqrtFlag = 0;
centerOnlyFlag = 0; % ONLY mean center input to PCA (dont divide by std)
centerNormalizeFlag = 1; % normalize cells/variables by their variance in addition to mean-centering (PCA on correlation matrix)    


allCellsPerTrialSpikeCount = data.allCellsPerTrialSpikeCount;
allTrials = data.allTrials;
bhvs2plot = data.bhvs2plot;
taxis = data.ax; 

% Single-Trial PCA 
% Input is nTrials x nCells , 
%    where each entry is # of spikes within a trial
%  Each point in state-space is a TRIAL 

if centerOnlyFlag
     [coeff,score,latent, ~] = pca(allCellsPerTrialSpikeCount); % centering (subtract column/cell means) done automatically by pca()
     neuralInput = allCellsPerTrialSpikeCount - mean(allCellsPerTrialSpikeCount); %for plotting only, later
elseif centerNormalizeFlag
     [coeff,score,latent, ~] = pca(allCellsPerTrialSpikeCount,'VariableWeights', 'variance'); 
     neuralInput = (allCellsPerTrialSpikeCount - mean(allCellsPerTrialSpikeCount)) ./ std(allCellsPerTrialSpikeCount,[],1);
elseif sqrtFlag 
    [coeff,score,latent, ~] = pca(sqrt(allCellsPerTrialSpikeCount)); 
    neuralInput = sqrt(allCellsPerTrialSpikeCount);
else
    [coeff,score,latent, ~] = pca(allCellsPerTrialSpikeCount,'Centered','off');
    neuralInput = allCellsPerTrialSpikeCount;
    disp('performing PCA on uncentered, non-normalized data!')
end


% explained variance 
explained = latent/sum(latent);

% eigenvalues 
figure; 
plot(cumsum(latent/sum(latent)), 'ok','MarkerFaceColor','k') 
xlabel('PC'); ylabel('% variance'); xlim([0 80]);
title('Cumulative Variance Explained');

% plot perTrial data projected onto first 3 PCs 
%%
[~,~,id] = unique(allTrials.trialType);
colors = 'rbgck';
markers = 'osd*.';
plotFig = figure;
for idx = 1:length(unique(allTrials.trialType))'
    dat2plot = score(allTrials.trialType == bhvs2plot(idx),:);
    scat(idx) = scatter3(dat2plot(:,1), dat2plot(:,2), dat2plot(:,3), 500, [colors(idx) markers(idx)], 'filled',...
        'MarkerFaceAlpha', 0.5,'MarkerEdgeColor',colors(idx),'linew',2);
    hold on;
end
xlabel('Neural PC1')
ylabel('Neural PC2')
zlabel('Neural PC3')
title(['Single Trials of Facial Gesture ']); 
legend('Threat','Lipsmack','Chew')

% % figure; imagesc(allCellsPerTrialSpikeCount'); colormap(jet); ylabel('cell'); xlabel('trial#'); 
% % colorbar;
% % xticks(1:length(allTrials.trialType));
% % xticklabels(allTrials.trialType);
% % yticks(1:size(allCellsPerTrialSpikeCount,2));
% % yticklabels(spikeLabels);
% % title('spikes PER trial, per Cell ')
% % 
% % 
% % figure; imagesc(neuralInput'); colormap(jet); ylabel('cell'); xlabel('trial#'); 
% % colorbar;
% % xticks(1:length(allTrials.trialType));
% % xticklabels(allTrials.trialType);
% % yticks(1:size(allCellsPerTrialSpikeCount,2));
% % yticklabels(spikeLabels);
% % title('spikes PER trial per cell, PREPROC INPUT ')
% % 
% % % % plot FR vs loadings to see if highFR neurons determine own PCs
% % figure;
% % subplot(1,2,1); plot(bin.meanFRs(bin.meanFRs >= threshlowFR), coeff(:,1), 'o'); grid; xlabel('meanFR'); ylabel('loadings on PC1')
% % [rho1, pval1] = corr(bin.meanFRs(bin.meanFRs >= threshlowFR)', coeff(:,1));
% % subplot(1,2,2); plot(bin.meanFRs(bin.meanFRs >= threshlowFR), coeff(:,2), 'o'); grid; xlabel('meanFR'); ylabel('loadings on PC2')
% % [rho2, pval2] = corr(bin.meanFRs(bin.meanFRs >= threshlowFR)', coeff(:,2));
% % sgtitle('Mean FRs vs. loadings on PC1, PC2')
% % 
% 
% % loadings of each neuron on each PC
%  figure; imagesc(coeff(:,1:5));  colorbar; colormap(bluewhitered)
%  ax = gca; 
%  ax.FontSize = 10;
%  yticks(1:size(allCellsPerTrialSpikeCount,2));
%  yticklabels(data.spikeLabels2plot);
%  xticks(1:5)
%  xlabel('PC','FontSize',18); ylabel('Cell ID','FontSize',18);
%  title('PC Loadings are Distributed Across Cells','FontSize',22)
% 

% % % biplots in 2D & 3D 
% % 
% % figure; h1 = biplot(coeff(:,1:3), 'scores', score(:,1:3), 'Color', 'b','varlabels', spikeLabels);
% % title(['biplot PC1-3; bhvtype: ' num2str(bhvs2plot)]); 
% % 
% % figure; h2 =  biplot(coeff(:,1:2), 'scores', score(:,1:2), 'Color', 'b','varlabels', spikeLabels);
% % title(['biplot PC1v2; bhvtype: ' num2str(bhvs2plot)]); 
% % 
% figure; h3 = biplot(coeff(:,[1,3]), 'scores', score(:,[1,3]), 'Color', 'b','varlabels', spikeLabels);
% title(['biplot PC1v3; bhvtype: ' num2str(bhvs2plot)]); 
% 
% % labeling neurons by array & trials by bhvType
% %colorBiplotVectorsByArray(allCellsPerTrialSpikeCount', spikeLabels2plot, h1, h2, h3);
% %colorObservationsOnBiplot(allCellsPerTrialSpikeCount', allTrials.trialType, h1, h2, h3)
% 