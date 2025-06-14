%%% this script with load the ANOVA on FR results 
%%% (which performed a 2-way anova (bhv, time, interaction) on FRs per cell
%%% This will calculate an omnibus statistic (chi-sq) to tell you if 
%%% 1) fraction of significantly modulated cells is same v. different
%%%     between regions
%%%
%%% 2) posthoc tests will tell you which region pairs had different
%%%     fractions of modulated cells 

clear all; 
subjects = {'barney','thor'};

pathToMatFiles = '/Users/geena/Dropbox/PhD/SUAinfo/';
pThres = 0.05;

%% get barney's data 
barneyRawData1 = load(fullfile(pathToMatFiles, 'Barney_210704/Data4Analysis/210704_barney_facePrefAnova.mat'));
barneyRawData2 = load(fullfile(pathToMatFiles, 'Barney_210706/Data4Analysis/210706_barney_facePrefAnova.mat'));
barneyRawData3 = load(fullfile(pathToMatFiles, 'Barney_210805/Data4Analysis/210805_barney_facePrefAnova.mat'));

barneyRawData.allResults = [barneyRawData1.allResults, barneyRawData2.allResults, barneyRawData3.allResults];
barneyRawData.allPvals = [barneyRawData1.allPvals, barneyRawData2.allPvals, barneyRawData3.allPvals];
barneyRawData.spikeLabels2plot = [barneyRawData1.spikeLabels2plot; barneyRawData2.spikeLabels2plot; barneyRawData3.spikeLabels2plot];

%% get thor's data 
thorRawData1 = load(fullfile(pathToMatFiles, 'Thor_171010/Data4Analysis/171010_Thor_facePrefAnova.mat'));
thorRawData2 = load(fullfile(pathToMatFiles, 'Thor_171027/Data4Analysis/171027_Thor_facePrefAnova.mat'));
thorRawData3 = load(fullfile(pathToMatFiles, 'Thor_171005/Data4Analysis/171005_Thor_facePrefAnova.mat'));
thorRawData4 = load(fullfile(pathToMatFiles, 'Thor_171128/Data4Analysis/171128_Thor_facePrefAnova.mat'));

thorRawData.allResults = [thorRawData1.allResults, thorRawData2.allResults, thorRawData3.allResults, thorRawData4.allResults];
thorRawData.allPvals = [thorRawData1.allPvals, thorRawData2.allPvals, thorRawData3.allPvals, thorRawData4.allPvals];
thorRawData.spikeLabels2plot = [thorRawData1.spikeLabels2plot; thorRawData2.spikeLabels2plot; thorRawData3.spikeLabels2plot; thorRawData4.spikeLabels2plot];

%% get the cortical regions per cell 
cortRegionOut = [];

for ss = 1:length(subjects)
    if strcmpi(subjects{ss},'Barney')
        data = barneyRawData;
    elseif strcmpi(subjects{ss},'Thor')
        data = thorRawData;
    end

    % pool all the p-values 
    for ii = 1:length(data.allResults)
        
        poolPvals(ii,:) = data.allPvals{1,ii};
        
    end
    poolPvalsOut{ss}(:,:) = poolPvals(~isnan(poolPvals(:,1)),:);
    spikeLabels2plot = data.spikeLabels2plot;
    
    % get region mapping
    [regions] = getChannel2CorticalRegionMapping(subjects{ss}, 1);
    
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
end

% threshold significance for two ANOVA factors + their interactions
% 
sigBhv = [poolPvalsOut{1}(:,1) ;poolPvalsOut{2}(:,1)] <= pThres;
sigTime = [poolPvalsOut{1}(:,2) ;poolPvalsOut{2}(:,2)] <= pThres;
sigInteract = [poolPvalsOut{1}(:,3) ;poolPvalsOut{2}(:,3)]  <= pThres;

% run the chi-square tests:
[tblBhv,chi2statBhv,pvalBhv] = crosstab(sigBhv,cortRegionOut);

[tblTime,chi2statTime,pvalTime] = crosstab(sigTime,cortRegionOut);

[tblInt,chi2statInt,pvalInt] = crosstab(sigInteract,cortRegionOut);

% now set-up the post-hoc tests 
% (adjusted pval is pThresh ./ number pairwise comparisons) 
pThresAdj = pThres ./ ((length(regions)*((length(regions))-1))/2);

S1Time = sigTime(cortRegionOut == 1)';
M1Time = sigTime(cortRegionOut == 2)';
PMvTime = sigTime(cortRegionOut == 3)';
M3Time = sigTime(cortRegionOut == 4)';

S1Int = sigInteract(cortRegionOut == 1)';
M1Int = sigInteract(cortRegionOut == 2)';
PMvInt = sigInteract(cortRegionOut == 3)';
M3Int = sigInteract(cortRegionOut == 4)';

% S1 vs M1 time
[tbl_S1M1time,chi2stat_S1M1time,pval_S1M1time,labels] = crosstab([S1Time M1Time],sort(cortRegionOut(cortRegionOut == 1 | cortRegionOut == 2)));

% S1 v PMv time
[tbl_S1PMvtime,chi2stat_S1PMvtime,pval_S1PMvtime] = crosstab([S1Time PMvTime],sort(cortRegionOut(cortRegionOut == 1 | cortRegionOut == 3)));

% S1 v M3 time
[tbl_S1M3time,chi2stat_S1M3time,pval_S1M3time] = crosstab([S1Time M3Time],sort(cortRegionOut(cortRegionOut == 1 | cortRegionOut == 4)));

% M1 v PMv time
[tbl_M1PMvtime,chi2stat_M1PMvtime,pval_M1PMvtime] = crosstab([M1Time PMvTime],sort(cortRegionOut(cortRegionOut == 2 | cortRegionOut == 3)));

% M1 v M3 time
[tbl_M1M3time,chi2stat_M1M3time,pval_M1M3time] = crosstab([M1Time M3Time],sort(cortRegionOut(cortRegionOut == 2 | cortRegionOut == 4)));

% PMv v M3 time
[tbl_PMvM3time,chi2stat_PMvM3time,pval_PMvM3time] = crosstab([PMvTime M3Time],sort(cortRegionOut(cortRegionOut == 3 | cortRegionOut == 4)));

%
%

% S1 vs M1 interact
[tbl_S1M1int,chi2stat_S1M1int,pval_S1M1int] = crosstab([S1Int M1Int],sort(cortRegionOut(cortRegionOut == 1 | cortRegionOut == 2)));

% S1 v PMv interact
[tbl_S1PMvint,chi2stat_S1PMvint,pval_S1PMvint] = crosstab([S1Int PMvInt],sort(cortRegionOut(cortRegionOut == 1 | cortRegionOut == 3)));

% S1 v M3 interact
[tbl_S1M3int,chi2stat_S1M3int,pval_S1M3int] = crosstab([S1Int M3Int],sort(cortRegionOut(cortRegionOut == 1 | cortRegionOut == 4)));

% M1 v PMv interact
[tbl_M1PMvint,chi2stat_M1PMvint,pval_M1PMvint] = crosstab([M1Int PMvInt],sort(cortRegionOut(cortRegionOut == 2 | cortRegionOut == 3)));

% M1 v M3 interact
[tbl_M1M3int,chi2stat_M1M3int,pval_M1M3int] = crosstab([M1Int M3Int],sort(cortRegionOut(cortRegionOut == 2 | cortRegionOut == 4)));

% PMv v M3 interact
[tbl_PMvM3int,chi2stat_PMvM3int,pval_PMvM3int] = crosstab([PMvInt M3Int],sort(cortRegionOut(cortRegionOut == 3 | cortRegionOut == 4)));




