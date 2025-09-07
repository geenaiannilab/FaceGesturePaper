% %% Simulate
rng(42);
% Blue bars = shuffle-corrected MI for each neuron.
% 
% Black error bars = 
% raw MI ± 1 SD of shuffle distribution.
% 
% Red crosses = 
% mean of the shuffle distribution (bias estimate).
clear all ; 

% dat = load('~/Desktop/binnedSpikes.mat');
% nTrials = groupcounts(dat.allTrials.trialType)'; 
% nGestures   = numel(nTrials);
% spikeCounts = squeeze(sum(dat.allSpikesAllSubs.binnedSpikes,2));
% nCells = size(spikeCounts,2);
% labels = dat.allTrials.trialType;

nTrials = [40, 100, 30];   
nGestures   = numel(nTrials);
nCells  = 2;

rates = [5  8 ;    
          10 9 ;   
          20 8 ];  
spikeCounts = [];
labels = [];

for c = 1:nGestures
    counts_c = poissrnd(repmat(rates(c,:), nTrials(c), 1));
    spikeCounts = [spikeCounts; counts_c];
    labels = [labels; c*ones(nTrials(c),1)];
end

%% Compute MI
nShuffles = 100;
[MI_corr, MI_raw, MI_shuff, pvals] = computeMI(spikeCounts, labels, nShuffles);

%% Display table
disp(table((1:nCells)', MI_raw', mean(MI_shuff)', MI_corr', pvals', ...
     'VariableNames', {'Neuron','MI_raw','ShuffleMean','MI_corr','pval'}));

%% Plot
figure; hold on;
bar(1:nCells, MI_corr, 'FaceColor',[0.2 0.6 0.8]); % corrected MI
errorbar(1:nCells, MI_raw, std(MI_shuff), 'k.', 'LineWidth',1.5); % raw ± shuffle SD
plot(1:nCells, mean(MI_shuff,1), 'rx', 'MarkerSize',10, 'LineWidth',2); % shuffle mean

% mark significant neurons (e.g. p<0.05) with a star
for c = 1:nCells
    if pvals(c) < 0.05
        text(c, MI_corr(c)+0.02, '*', 'HorizontalAlignment','center','FontSize',14);
    end
end

xlabel('Neuron');
ylabel('Mutual Information (bits)');
legend({'Bias-corrected MI','Raw ± shuffle SD','Shuffle mean'}, 'Location','NorthWest');
title('Per-cell MI with Shuffle Correction and Significance');
set(gca,'FontSize',12);
