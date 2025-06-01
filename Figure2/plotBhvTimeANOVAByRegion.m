clear all
set(0,'defaultAxesFontSize',20)

%% Control options
includeF4 = true;


%% Load matfiles
pathToMatFiles = '/Users/geena/Dropbox/PhD/SUAinfo';
thorRawData1 = load(fullfile(pathToMatFiles, 'Thor_171010/Data4Analysis/171010_Thor_facePrefAnova.mat'));
thorRawData2 = load(fullfile(pathToMatFiles, 'Thor_171027/Data4Analysis/171027_Thor_facePrefAnova.mat'));
thorRawData3 = load(fullfile(pathToMatFiles, 'Thor_171005/Data4Analysis/171005_Thor_facePrefAnova.mat'));
thorRawData4 = load(fullfile(pathToMatFiles, 'Thor_171128/Data4Analysis/171128_Thor_facePrefAnova.mat'));

thorRawData.allResults = [thorRawData1.allResults, thorRawData2.allResults, thorRawData3.allResults, thorRawData4.allResults];
thorRawData.allPvals = [thorRawData1.allPvals, thorRawData2.allPvals, thorRawData3.allPvals, thorRawData4.allPvals];
thorRawData.spikeLabels2plot = [thorRawData1.spikeLabels2plot; thorRawData2.spikeLabels2plot; thorRawData3.spikeLabels2plot; thorRawData4.spikeLabels2plot];


barneyRawData1 = load(fullfile(pathToMatFiles, 'Barney_210704/Data4Analysis/210704_barney_facePrefAnova.mat'));
barneyRawData2 = load(fullfile(pathToMatFiles, 'Barney_210706/Data4Analysis/210706_barney_facePrefAnova.mat'));
barneyRawData3 = load(fullfile(pathToMatFiles, 'Barney_210805/Data4Analysis/210805_barney_facePrefAnova.mat'));

barneyRawData.allResults = [barneyRawData1.allResults, barneyRawData2.allResults,barneyRawData3.allResults ];
barneyRawData.allPvals = [barneyRawData1.allPvals, barneyRawData2.allPvals, barneyRawData3.allPvals];
barneyRawData.spikeLabels2plot = [barneyRawData1.spikeLabels2plot; barneyRawData2.spikeLabels2plot; barneyRawData3.spikeLabels2plot];

%% Set up basic variables
pValueThreshold = 0.05;


%% Get Barney's data

% Remove NaNs
for ii = 1:length(barneyRawData.allResults)
    
    poolPvals(ii,:) = barneyRawData.allPvals{1,ii};
    
end
poolPvals = poolPvals(~isnan(poolPvals(:,1)),:);
spikeLabels2plot = barneyRawData.spikeLabels2plot;

% Define spikeLabels2plot
regions{1}.label = 'S1';
regions{1}.channels = 1:32;
regions{2}.label = 'M1';
regions{2}.channels = 33:128;
% regions{2}.label = 'M1m';
% regions{2}.channels = 33:64;
% if includeF4
%     regions{3}.label = 'M1lat';
%     regions{3}.channels = 65:128;
% else
%     regions{3}.label = 'M1lat';
%     regions{3}.channels = 65:96;
% end
regions{3}.label = 'PMv';
regions{3}.channels = 129:192;
regions{4}.label = 'M3';
regions{4}.channels = 193:240;

factors = {'behavior', 'time', 'interaction'};

% Loop around regions. For each region: 1) identify relevant indices from
% spikeLabels2Plot, 2) take these indices and count the number of
% corresponding cells with a significant p-value, and 3) stash the result
for ff = 1:length(factors)
    for rr = 1:length(regions)
        regionIndices = find(spikeLabels2plot(:,1) >= min(regions{rr}.channels) & spikeLabels2plot(:,1) <= max(regions{rr}.channels));
        totalNumberOfCells = length(regionIndices);
        significantCells = find(poolPvals(regionIndices,ff) <= pValueThreshold);
        
        barneyResults.(factors{ff}).(regions{rr}.label).percentSignificant = 100*(length(significantCells)/totalNumberOfCells);
        barneyResults.(factors{ff}).(regions{rr}.label).numberSignificant = length(significantCells);
        barneyResults.(factors{ff}).(regions{rr}.label).totalNumberOfCells = totalNumberOfCells;
        
    end
    
    %barneyResults.(factors{ff}).lateralRegions.percentSignificant = (barneyResults.(factors{ff}).M1.percentSignificant * barneyResults.(factors{ff}).M1.totalNumberOfCells +  barneyResults.(factors{ff}).PMv.percentSignificant * barneyResults.(factors{ff}).PMv.totalNumberOfCells)/...
   %     (barneyResults.(factors{ff}).M1.totalNumberOfCells + barneyResults.(factors{ff}).PMv.totalNumberOfCells);
    
%     barneyResults.(factors{ff}).combinedM1.percentSignificant = (barneyResults.(factors{ff}).M1.percentSignificant * barneyResults.(factors{ff}).M1.totalNumberOfCells + barneyResults.(factors{ff}).M1lat.percentSignificant * barneyResults.(factors{ff}).M1lat.totalNumberOfCells)/...
%         (barneyResults.(factors{ff}).M1m.totalNumberOfCells + barneyResults.(factors{ff}).M1lat.totalNumberOfCells);
%     
 %   barneyResults.(factors{ff}).nonM1.percentSignificant = (barneyResults.(factors{ff}).M3.percentSignificant * barneyResults.(factors{ff}).M3.totalNumberOfCells + barneyResults.(factors{ff}).PMv.percentSignificant * barneyResults.(factors{ff}).PMv.totalNumberOfCells)/...
%        (barneyResults.(factors{ff}).M3.totalNumberOfCells + barneyResults.(factors{ff}).PMv.totalNumberOfCells);

    
end

%% Get Thor's data


% Remove NaNs
for ii = 1:length(thorRawData.allResults)
    
    poolPvals(ii,:) = thorRawData.allPvals{1,ii};
    
end
poolPvals = poolPvals(~isnan(poolPvals(:,1)),:);
spikeLabels2plot = thorRawData.spikeLabels2plot;

% Define spikeLabels2plot
regions{1}.label = 'S1';
regions{1}.channels = 1:32;
regions{2}.label = 'M1';
regions{2}.channels = 33:96;
%regions{2}.label = 'M1m';
%regions{2}.channels = 33:64;
%regions{3}.label = 'M1lat';
%regions{3}.channels = 65:96;
regions{3}.label = 'PMv';
regions{3}.channels = 97:128;
regions{4}.label = 'M3';
regions{4}.channels = 129:192;

factors = {'behavior', 'time', 'interaction'};

% Loop around regions. For each region: 1) identify relevant indices from
% spikeLabels2Plot, 2) take these indices and count the number of
% corresponding cells with a significant p-value, and 3) stash the result
for ff = 1:length(factors)
    for rr = 1:length(regions)
        regionIndices = find(spikeLabels2plot(:,1) >= min(regions{rr}.channels) & spikeLabels2plot(:,1) <= max(regions{rr}.channels));
        totalNumberOfCells = length(regionIndices);
        significantCells = find(poolPvals(regionIndices,ff) <= pValueThreshold);
        
        thorResults.(factors{ff}).(regions{rr}.label).percentSignificant = 100*(length(significantCells)/totalNumberOfCells);
        thorResults.(factors{ff}).(regions{rr}.label).numberSignificant = length(significantCells);
        thorResults.(factors{ff}).(regions{rr}.label).totalNumberOfCells = totalNumberOfCells;
        
    end
    %thorResults.(factors{ff}).lateralRegions.percentSignificant = (thorResults.(factors{ff}).M1m.percentSignificant * thorResults.(factors{ff}).M1m.totalNumberOfCells + thorResults.(factors{ff}).M1lat.percentSignificant * thorResults.(factors{ff}).M1lat.totalNumberOfCells)/...
     %   (thorResults.(factors{ff}).M1m.totalNumberOfCells + thorResults.(factors{ff}).M1lat.totalNumberOfCells);
    %thorResults.(factors{ff}).combinedM1.percentSignificant = (thorResults.(factors{ff}).M1m.percentSignificant * thorResults.(factors{ff}).M1m.totalNumberOfCells + thorResults.(factors{ff}).M1lat.percentSignificant * thorResults.(factors{ff}).M1lat.totalNumberOfCells)/...
    %    (thorResults.(factors{ff}).M1m.totalNumberOfCells + thorResults.(factors{ff}).M1lat.totalNumberOfCells);
    %   thorResults.(factors{ff}).nonM1.percentSignificant = (thorResults.(factors{ff}).M3.percentSignificant * thorResults.(factors{ff}).M3.totalNumberOfCells + thorResults.(factors{ff}).PMv.percentSignificant * thorResults.(factors{ff}).PMv.totalNumberOfCells)/...
    %    (thorResults.(factors{ff}).M3.totalNumberOfCells + thorResults.(factors{ff}).PMv.totalNumberOfCells);

end

for ff = 1:length(factors)
    for rr = 1:length(regions)
        combinedResults.(factors{ff}).(regions{rr}.label).numberSignificant = thorResults.(factors{ff}).(regions{rr}.label).numberSignificant + barneyResults.(factors{ff}).(regions{rr}.label).numberSignificant;
        combinedResults.(factors{ff}).(regions{rr}.label).totalNumberOfCells  = thorResults.(factors{ff}).(regions{rr}.label).totalNumberOfCells + barneyResults.(factors{ff}).(regions{rr}.label).totalNumberOfCells;
        combinedResults.(factors{ff}).(regions{rr}.label).percentSignificant = 100*(combinedResults.(factors{ff}).(regions{rr}.label).numberSignificant / combinedResults.(factors{ff}).(regions{rr}.label).totalNumberOfCells);

    end
end

%% Plot
plotFig = figure;
for ff = 1:length(factors)
    subplot(1,3,ff)
    b = bar([barneyResults.(factors{ff}).(regions{1}.label).percentSignificant, barneyResults.(factors{ff}).(regions{2}.label).percentSignificant, barneyResults.(factors{ff}).(regions{3}.label).percentSignificant, barneyResults.(factors{ff}).(regions{4}.label).percentSignificant; ...
        thorResults.(factors{ff}).(regions{1}.label).percentSignificant, thorResults.(factors{ff}).(regions{2}.label).percentSignificant, thorResults.(factors{ff}).(regions{3}.label).percentSignificant, thorResults.(factors{ff}).(regions{4}.label).percentSignificant]', 'FaceColor','flat');
    b(1).CData = [         
    0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	
         ];
    b(2).CData =  [         0.4940 0.1840 0.5560
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	
        ];
   
    b(2).FaceAlpha = 0.5;

    title(factors{ff})
    xticklabels({regions{1}.label, regions{2}.label, regions{3}.label, regions{4}.label})
    xlabel('Region')
    ylim([0 90]);
    ylabel('% Significant Cells')
%     if ff == 3
%         legend('Monkey B', 'Monkey T')
%     end
end




%% Plot
plotFig = figure;
for ff = 1:length(factors)
    subplot(1,3,ff)
    b = bar([combinedResults.(factors{ff}).(regions{1}.label).percentSignificant, combinedResults.(factors{ff}).(regions{2}.label).percentSignificant, combinedResults.(factors{ff}).(regions{3}.label).percentSignificant, combinedResults.(factors{ff}).(regions{4}.label).percentSignificant]', 'FaceColor','flat');
    b(1).CData = [         
    0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	
         ];

    title(factors{ff})
    xticklabels({regions{1}.label, regions{2}.label, regions{3}.label, regions{4}.label})
    xlabel('Region')
    ylim([0 90]);
    ylabel('% Significant Cells')
    
end



plotFig = figure;
for ff = 1:length(factors)
    subplot(1,3,ff)
    b = bar([barneyResults.(factors{ff}).lateralRegions.percentSignificant, barneyResults.(factors{ff}).M3.percentSignificant; ...
        thorResults.(factors{ff}).lateralRegions.percentSignificant, thorResults.(factors{ff}).M3.percentSignificant]');
    title(factors{ff})
    xticklabels({'Lateral', 'Medial'})
    xlabel('Region')
    ylim([0 1]);
    ylabel('% Significant Cells')
    if ff == 3
        legend('Monkey B', 'Monkey T')
    end
end

plotFig = figure;
for ff = 1:length(factors)
    subplot(1,3,ff)
    b = bar([barneyResults.(factors{ff}).combinedM1.percentSignificant, barneyResults.(factors{ff}).nonM1.percentSignificant; ...
        thorResults.(factors{ff}).combinedM1.percentSignificant, thorResults.(factors{ff}).nonM1.percentSignificant]');
    title(factors{ff})
    xticklabels({'M1', 'M3+PMv'})
    xlabel('Region')
    ylim([0 1]);
    ylabel('% Significant Cells')
    if ff == 3
        legend('Monkey B', 'Monkey T')
    end
end