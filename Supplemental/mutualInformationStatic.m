% %% Simulate
% %rng(42);
% Blue bars = shuffle-corrected MI for each neuron.
% 
% Black error bars = 
% raw MI ± 1 SD of shuffle distribution.
% 
% Red crosses = 
% mean of the shuffle distribution (bias estimate).
clear all ; close all 

%% Setup
dates = {'Barney_210704','Barney_210706','Barney_210805', ...
         'Thor_171005','Thor_171010','Thor_171027','Thor_171128'};
baseDir = '~/Dropbox/PhD/SUAInfo/';
saveFlag = 1;

binSize   = 0.001;  % sec
nShuffles = 10000;
alpha     = 0.05;
time2keep = [-1 1];

% for plotting
colorMap = [         
    0.4940 0.1840 0.5560
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980
    0.9290 0.6940 0.1250 ];

% Pooled results across days
MI_corr_all = [];
MI_raw_all  = [];
pvals_all   = [];
regions_all = {};
day_all     = {};

% Per-day tidy table (grow as we go)
PerDay = table('Size',[0 9], ...
    'VariableTypes', {'string','string','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Day','Region','nCells','nSig','FracSig','MIcorr_mean','MIcorr_median','MIraw_mean','MIraw_median'});

%% Main loop over days
for dd = 1:numel(dates)
    dayID = dates{dd};
    workDir = [baseDir dayID '/Data4Analysis/'];
    dat = load([workDir 'binnedSpikesForMI.mat']); 
    originalTimeBase = dat.taxis2take;
    tIndx2keep = (originalTimeBase >= time2keep(1) & originalTimeBase <= time2keep(2));
    
    % --- get day-specific data ---
    spikes = dat.allSpikesAllSubs.binnedSpikes(:,tIndx2keep,:); % [nTrials x nTime x nCells]
    labels = dat.allTrials.trialType;

    % --- preprocess: remove near-silent cells ---
    cellTotals = squeeze(sum(sum(spikes,1),2));
    keepCells  = cellTotals > 1;
    spikes = spikes(:,:,keepCells);
    regions_day = dat.cortRegion(keepCells);

    % collapse across time bins for MI (non-time-resolved)
    spikeCounts = squeeze(sum(spikes,2)); % [nTrials x nCells_day]

    % --- compute MI ---
    [MI_corr, MI_raw, MI_shuff, pvals] = computeMI(spikeCounts, labels, nShuffles);

    % pool (per-cell)
    MI_corr_all = [MI_corr_all MI_corr];
    MI_raw_all  = [MI_raw_all  MI_raw];
    pvals_all   = [pvals_all   pvals];
    regions_all = [regions_all regions_day(:)'];
    day_all     = [day_all repmat({dayID},1,numel(regions_day))];

    % ---- per-day × region summary rows ----
    uR = unique(regions_day,'stable');
    isSig_day = pvals < alpha;
    for r = 1:numel(uR)
        idx = strcmp(regions_day, uR{r});
        if ~any(idx), continue; end
        nTot = sum(idx);
        nSig = sum(isSig_day(idx));
        PerDay = [PerDay; {
            string(dayID), string(uR{r}), ...
            nTot, nSig, nSig/max(nTot,1), ...
            mean(MI_corr(idx)), median(MI_corr(idx)), ...
            mean(MI_raw(idx)),  median(MI_raw(idx))
        }];
    end
end

%% Region-wise pooled summaries & stats
regions = regions_all;
uRegs   = unique(regions,'stable');
nRegs   = numel(uRegs);

% fraction significant per region (pooled)
isSig   = (pvals_all < alpha);
fracSig = zeros(nRegs,1);
counts  = zeros(nRegs,2); % [#sig, #nonsig]
for i = 1:nRegs
    idx = strcmp(regions, uRegs{i});
    nSig = sum(isSig(idx));
    nTot = sum(idx);
    fracSig(i) = nSig / max(nTot,1);
    counts(i,:) = [nSig, max(nTot-nSig,0)];
end

% --- Chi-square test of independence on pooled fractions ---
[chi2stat,p_chi,stats_chi] = chi2indep(counts);
fprintf('Chi-square (pooled) across regions: χ²(%d) = %.2f, p = %.4f\n', ...
        stats_chi.df, chi2stat, p_chi);

%%%
if saveFlag 
    save([baseDir 'MIResults_allDays.mat']);
end
%%%%%%%%%%%%%%%%%%%%%


%% --- Helper function to overlay sig bars ---
function overlaySig(results)
    ylim_curr = ylim;
    y_top = ylim_curr(2);
    step = 0.05 * range(ylim_curr); % vertical spacing
    for r = 1:size(results,1)
        g1 = results(r,1); g2 = results(r,2); pval = results(r,6);
        if pval < 0.05
            addSigBar(g1,g2,y_top + r*step,pval);
        end
    end
    ylim([ylim_curr(1) y_top + (size(results,1)+2)*step]);
end

%% --- Draw a bar with stars ---
function addSigBar(x1,x2,y,p)
    line([x1 x2],[y y],'Color','k','LineWidth',1.5);
    line([x1 x1],[y-0.01 y],'Color','k','LineWidth',1.5);
    line([x2 x2],[y-0.01 y],'Color','k','LineWidth',1.5);
    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    else
        stars = '*';
    end
    text(mean([x1 x2]), y+0.01, stars, ...
        'HorizontalAlignment','center','FontSize',12);
end

function [chi2stat, p, stats] = chi2indep(counts)
    % counts: contingency table (rows = groups, cols = categories)
    % returns chi-square statistic, p-value, and stats struct

    rowSums = sum(counts,2);
    colSums = sum(counts,1);
    N       = sum(rowSums);

    expected = (rowSums * colSums) / N;
    chi2stat = sum((counts - expected).^2 ./ expected,'all');

    df = (size(counts,1)-1) * (size(counts,2)-1);
    p  = 1 - chi2cdf(chi2stat, df);

    stats.df = df;
    stats.expected = expected;
end
