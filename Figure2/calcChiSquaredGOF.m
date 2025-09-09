%%% this script with load the ANOVA on FR results 
%%% (which performed a 2-way anova (bhv, time, interaction) on FRs per cell
%%% This will calculate an omnibus statistic (chi-sq) to tell you if 
%%% 1) fraction of significantly modulated cells is same v. different
%%%     between regions
%%%
%%% 2) posthoc tests will tell you which region pairs had different
%%%     fractions of modulated cells 

function [results] = calcChiSquaredGOF(data,subj)

% thorRawData.allResults = [thorRawData1.allResults, thorRawData2.allResults, thorRawData3.allResults, thorRawData4.allResults];
% thorRawData.allPvals = [thorRawData1.allPvals, thorRawData2.allPvals, thorRawData3.allPvals, thorRawData4.allPvals];
% thorRawData.spikeLabels2plot = [thorRawData1.spikeLabels2plot; thorRawData2.spikeLabels2plot; thorRawData3.spikeLabels2plot; thorRawData4.spikeLabels2plot];

%% get the cortical regions per cell 
cortRegionOut = [];

% pool all the p-values
for ii = 1:length(data.allResults)

    poolPvals(ii,:) = data.allPvals{1,ii};

end
poolPvalsOut{ss}(:,:) = poolPvals(~isnan(poolPvals(:,1)),:);
spikeLabels2plot = data.spikeLabels2plot;

% get region mapping
[regions] = getChannel2CorticalRegionMapping(subj, 1);

% code the cortical region for each cell
for unit = 1:length(spikeLabels2plot)
    el = spikeLabels2plot(unit,1);

    if ismember(el, regions{1,1}.channels)
        cortRegion(unit) = 1; % S1
    elseif ismember(el, regions{1,2}.channels)
        cortRegion(unit) = 2; % M1
    elseif ismember(el, regions{1,3}.channels)
        cortRegion(unit) = 3; % PMv
    elseif ismember(el, regions{1,4}.channels)
        cortRegion(unit) = 4;
    else
        disp('warning!');
    end

end

cortRegionOut = [cortRegionOut; cortRegion'];

clear cortRegion;
clear poolPvals;

% threshold significance for two ANOVA factors + their interactions
%
sigBhv = [poolPvalsOut{1}(:,1) ;poolPvalsOut{2}(:,1)] <= pThres;
sigTime = [poolPvalsOut{1}(:,2) ;poolPvalsOut{2}(:,2)] <= pThres;
sigInteract = [poolPvalsOut{1}(:,3) ;poolPvalsOut{2}(:,3)]  <= pThres;

% run the chi-square tests:
[tblBhv,chi2statBhv,pvalBhv] = crosstab(sigBhv,cortRegionOut);

[tblTime,chi2statTime,pvalTime] = crosstab(sigTime,cortRegionOut);

[tblInt,chi2statInt,pvalInt] = crosstab(sigInteract,cortRegionOut);


results.bhv.table = tblBhv;
results.bhv.chi2stat = chi2statBhv;
results.bhv.pval = pvalBhv;

results.time.table = tblTime;
results.time.chi2stat = chi2statTime;
results.time.pval = pvalTime;

results.int.table = tblInt;
results.int.chi2stat = chi2statInt;
results.int.pval = pvalInt;

%% ---- Package pairwise tests into results struct ----
% Meta
regionNames = {'S1','M1','PMv','M3'};
nPairs      = numel(regionNames)*(numel(regionNames)-1)/2;
if ~exist('pThres','var'), pThres = 0.05; end
if ~exist('pThresAdj','var')
    pThresAdj = pThres / nPairs;  % Bonferroni familywise alpha for pairwise tests
end

% Helper to build an entry
mk = @(tbl,chi2,p,r1,r2) struct( ...
    'regions',      {{r1,r2}}, ...
    'table',        tbl, ...
    'chi2stat',     chi2, ...
    'pval',         p, ...
    'pval_bonf',    min(p*nPairs,1), ...   % Bonferroni-adjusted p
    'alpha',        pThres, ...
    'alpha_bonf',   pThresAdj, ...
    'sig_raw',      p < pThres, ...
    'sig_bonf',     p < pThresAdj );

% TIME pairwise
results.time.alpha       = pThres;
results.time.alpha_bonf  = pThresAdj;
results.time.pairwise = struct();
results.time.pairwise.S1_vs_M1 = mk(tbl_S1M1time,  chi2stat_S1M1time,  pval_S1M1time,  'S1','M1');
results.time.pairwise.S1_vs_PMv= mk(tbl_S1PMvtime, chi2stat_S1PMvtime, pval_S1PMvtime, 'S1','PMv');
results.time.pairwise.S1_vs_M3 = mk(tbl_S1M3time,  chi2stat_S1M3time,  pval_S1M3time,  'S1','M3');
results.time.pairwise.M1_vs_PMv= mk(tbl_M1PMvtime, chi2stat_M1PMvtime, pval_M1PMvtime, 'M1','PMv');
results.time.pairwise.M1_vs_M3 = mk(tbl_M1M3time,  chi2stat_M1M3time,  pval_M1M3time,  'M1','M3');
results.time.pairwise.PMv_vs_M3= mk(tbl_PMvM3time, chi2stat_PMvM3time, pval_PMvM3time, 'PMv','M3');

% INTERACTION pairwise
results.int.alpha       = pThres;
results.int.alpha_bonf  = pThresAdj;
results.int.pairwise = struct();
results.int.pairwise.S1_vs_M1 = mk(tbl_S1M1int,  chi2stat_S1M1int,  pval_S1M1int,  'S1','M1');
results.int.pairwise.S1_vs_PMv= mk(tbl_S1PMvint, chi2stat_S1PMvint, pval_S1PMvint, 'S1','PMv');
results.int.pairwise.S1_vs_M3 = mk(tbl_S1M3int,  chi2stat_S1M3int,  pval_S1M3int,  'S1','M3');
results.int.pairwise.M1_vs_PMv= mk(tbl_M1PMvint, chi2stat_M1PMvint, pval_M1PMvint, 'M1','PMv');
results.int.pairwise.M1_vs_M3 = mk(tbl_M1M3int,  chi2stat_M1M3int,  pval_M1M3int,  'M1','M3');
results.int.pairwise.PMv_vs_M3= mk(tbl_PMvM3int, chi2stat_PMvM3int, pval_PMvM3int, 'PMv','M3');

% (Optional) also keep region list & nPairs for provenance
results.meta.regionNames = regionNames;
results.meta.nPairs      = nPairs;
end