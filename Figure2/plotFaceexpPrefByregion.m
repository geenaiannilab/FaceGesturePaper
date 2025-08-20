clear all
set(0,'defaultAxesFontSize',20)

%% Control options
includeF4 = true;


%% Load matfiles
pathToMatFiles = '/Users/geena/Dropbox/PhD/SUAinfo';
thorRawData1 = load(fullfile(pathToMatFiles, 'Thor_171010/Data4Analysis/171010_Thor_faceExpPref_V2.mat'));
thorRawData2 = load(fullfile(pathToMatFiles, 'Thor_171027/Data4Analysis/171027_Thor_faceExpPref_V2.mat'));
thorRawData3 = load(fullfile(pathToMatFiles, 'Thor_171005/Data4Analysis/171005_Thor_faceExpPref_V2.mat'));
thorRawData4 = load(fullfile(pathToMatFiles, 'Thor_171128/Data4Analysis/171128_Thor_faceExpPref_V2.mat'));

thorRawData.depthPref = [thorRawData1.depthPref, thorRawData2.depthPref, thorRawData3.depthPref,thorRawData4.depthPref];
thorRawData.faceExpPref = [thorRawData1.faceExpPref, thorRawData2.faceExpPref, thorRawData3.faceExpPref,thorRawData4.faceExpPref];
thorRawData.spikeLabels2plot = [thorRawData1.spikeLabels2plot; thorRawData2.spikeLabels2plot; thorRawData3.spikeLabels2plot; thorRawData4.spikeLabels2plot];


barneyRawData1 = load(fullfile(pathToMatFiles, 'Barney_210704/Data4Analysis/210704_barney_faceExpPref_V2.mat'));
barneyRawData2 = load(fullfile(pathToMatFiles, 'Barney_210805/Data4Analysis/210805_barney_faceExpPref_V2.mat'));
barneyRawData3 = load(fullfile(pathToMatFiles, 'Barney_210706/Data4Analysis/210706_barney_faceExpPref_V2.mat'));

barneyRawData.depthPref = [barneyRawData1.depthPref, barneyRawData2.depthPref, barneyRawData3.depthPref];
barneyRawData.faceExpPref = [barneyRawData1.faceExpPref, barneyRawData2.faceExpPref, barneyRawData3.faceExpPref];
barneyRawData.spikeLabels2plot = [barneyRawData1.spikeLabels2plot; barneyRawData2.spikeLabels2plot; barneyRawData3.spikeLabels2plot];

%% Set up basic variables
pValueThreshold = 0.05;


%% Get Barney's data

spikeLabels2plot = barneyRawData.spikeLabels2plot;

% Define spikeLabels2plot
regions{1}.label = 'S1';
regions{1}.channels = 1:32;
regions{2}.label = 'M1';
regions{2}.channels = 33:128;
regions{3}.label = 'PMv';
regions{3}.channels = 129:192;
regions{4}.label = 'M3';
regions{4}.channels = 193:240;

% regions{2}.label = 'M1m';
% regions{2}.channels = 33:64;
% if includeF4
%     regions{3}.label = 'M1lat';
%     regions{3}.channels = 65:128;
% else
%     regions{3}.label = 'M1lat';
%     regions{3}.channels = 65:96;
% end
% regions{4}.label = 'PMv';
% regions{4}.channels = 129:192;
% regions{5}.label = 'M3';
% regions{5}.channels = 193:240;


% Loop around regions. For each region: 1) identify relevant indices from
% spikeLabels2Plot, 2) take these indices and & stash the result
    for rr = 1:length(regions)
        regionIndices = find(spikeLabels2plot(:,1) >= min(regions{rr}.channels) & spikeLabels2plot(:,1) <= max(regions{rr}.channels));
        totalNumberOfCells = length(regionIndices);
        barneyResults.(regions{rr}.label).faceExpPrefIndex = barneyRawData.faceExpPref(regionIndices);
        barneyResults.(regions{rr}.label).depthPrefIndex = barneyRawData.depthPref(regionIndices);
        barneyResults.(regions{rr}.label).totalNumberOfCells = totalNumberOfCells;
        
    end
    
% put together all the results 
% barneyResults.lateralRegions.faceExpPrefIndex = cat(2, barneyResults.M1m.faceExpPrefIndex, barneyResults.M1lat.faceExpPrefIndex, barneyResults.PMv.faceExpPrefIndex);
% barneyResults.lateralRegions.depthPrefIndex = cat(2, barneyResults.M1m.depthPrefIndex, barneyResults.M1lat.depthPrefIndex, barneyResults.PMv.depthPrefIndex);
% barneyResults.lateralRegions.totalNumberOfCells = barneyResults.M1m.totalNumberOfCells + barneyResults.M1lat.totalNumberOfCells + barneyResults.PMv.totalNumberOfCells;
% 
barneyResults.combinedM1.faceExpPrefIndex = cat(2, barneyResults.M1.faceExpPrefIndex);
barneyResults.combinedM1.depthPrefIndex = cat(2, barneyResults.M1.depthPrefIndex);
barneyResults.combinedM1.totalNumberOfCells = barneyResults.M1.totalNumberOfCells;

barneyResults.nonM1.faceExpPrefIndex = cat(2, barneyResults.PMv.faceExpPrefIndex, barneyResults.M3.faceExpPrefIndex);
barneyResults.nonM1.depthPrefIndex = cat(2, barneyResults.PMv.depthPrefIndex, barneyResults.M3.depthPrefIndex);
barneyResults.nonM1.totalNumberOfCells = barneyResults.PMv.totalNumberOfCells + barneyResults.M3.totalNumberOfCells;



%% Get Thor's data

spikeLabels2plot = thorRawData.spikeLabels2plot;

% Define spikeLabels2plot
regions{1}.label = 'S1';
regions{1}.channels = 1:32;
regions{2}.label = 'M1';
regions{2}.channels = 33:96;
regions{3}.label = 'PMv';
regions{3}.channels = 97:128;
regions{4}.label = 'M3';
regions{4}.channels = 129:192;

% regions{2}.label = 'M1m';
% regions{2}.channels = 33:64;
% regions{3}.label = 'M1lat';
% regions{3}.channels = 65:96;
% regions{4}.label = 'PMv';
% regions{4}.channels = 97:128;
% regions{5}.label = 'M3';
% regions{5}.channels = 129:192;


% Loop around regions. For each region: 1) identify relevant indices from
% spikeLabels2Plot, & stash the result
for rr = 1:length(regions)
    regionIndices = find(spikeLabels2plot(:,1) >= min(regions{rr}.channels) & spikeLabels2plot(:,1) <= max(regions{rr}.channels));
    totalNumberOfCells = length(regionIndices);
    thorResults.(regions{rr}.label).faceExpPrefIndex = thorRawData.faceExpPref(regionIndices);
    thorResults.(regions{rr}.label).depthPrefIndex = thorRawData.depthPref(regionIndices);
    thorResults.(regions{rr}.label).totalNumberOfCells = totalNumberOfCells;
end
  
for rr = 1:length(regions)
    
    combinedResults.(regions{rr}.label).faceExpPrefIndex = [barneyResults.(regions{rr}.label).faceExpPrefIndex thorResults.(regions{rr}.label).faceExpPrefIndex];
    combinedResults.(regions{rr}.label).depthPrefIndx = [barneyResults.(regions{rr}.label).depthPrefIndex thorResults.(regions{rr}.label).depthPrefIndex];
    combinedResults.(regions{rr}.label).totalNumberOfCells = barneyResults.(regions{rr}.label).totalNumberOfCells + thorResults.(regions{rr}.label).totalNumberOfCells;
end

% put together all the results 
% thorResults.lateralRegions.faceExpPrefIndex = cat(2, thorResults.M1m.faceExpPrefIndex, thorResults.M1lat.faceExpPrefIndex, thorResults.PMv.faceExpPrefIndex);
% thorResults.lateralRegions.depthPrefIndex = cat(2, thorResults.M1m.depthPrefIndex, thorResults.M1lat.depthPrefIndex, thorResults.PMv.depthPrefIndex);
% thorResults.lateralRegions.totalNumberOfCells = thorResults.M1m.totalNumberOfCells + thorResults.M1lat.totalNumberOfCells + thorResults.PMv.totalNumberOfCells;

thorResults.combinedM1.faceExpPrefIndex = cat(2, thorResults.M1.faceExpPrefIndex);
thorResults.combinedM1.depthPrefIndex = cat(2, thorResults.M1.depthPrefIndex);
thorResults.combinedM1.totalNumberOfCells = thorResults.M1.totalNumberOfCells ;

thorResults.nonM1.faceExpPrefIndex = cat(2, thorResults.PMv.faceExpPrefIndex, thorResults.M3.faceExpPrefIndex);
thorResults.nonM1.depthPrefIndex = cat(2, thorResults.PMv.depthPrefIndex, thorResults.M3.depthPrefIndex);
thorResults.nonM1.totalNumberOfCells = thorResults.PMv.totalNumberOfCells + thorResults.M3.totalNumberOfCells;

% concatente both animals for plotting 
cmap1 = [         
    0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	
         ];
cmap2 = [         0.4940 0.1840 0.5560
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	
        ];
    b(2).FaceAlpha = 0.5;

 xall = [];
 yall = [];
 colorAll = [];
 faceAlphaAll = [];
 counter = 1;
for ii = 1:length(regions)
    for subj = 1:2
        if subj == 1
            thisRegion = barneyResults.(regions{ii}.label);
            color = cmap1(ii,:);
            faceAlpha = 0.8;
            %color = [0 0.4470 0.7410];
        elseif subj == 2
            thisRegion = thorResults.(regions{ii}.label);
            color = cmap1(ii,:);
            faceAlpha = 0.5;
            %color = [0.8500 0.3250 0.0980];
        end
        x = ones(thisRegion.totalNumberOfCells,1)*(counter);
        y = thisRegion.faceExpPrefIndex';

        xall = cat(1,xall, x);
        yall = cat(1,yall, y);
        colorAll = cat(1,colorAll, color);
        faceAlphaAll = cat(1,faceAlphaAll, faceAlpha);
        if subj == 1
            counter = counter +1;
        elseif subj == 2
            counter = counter +1.75;
        end
    end
end

%% Plot 
plotFig = figure;
b = beeswarm(xall,yall,'Colormap', colorAll,'dot_size',1.2,'MarkerFaceAlpha',faceAlphaAll','overlay_style','ci','sort_style','fan');
xticks([1.5 4.25 7 9.75 ]);
xticklabels({'S1', 'M1','PMv', 'M3'})
xlabel('Region')
ylim([0 1.1]);
ylabel('Preference Index');
title('Face Expression Preference Indices per Regions')
legend('Monkey B','','','Monkey T');

%% both animals together 
% concatente both animals for plotting 
cmap1 = [         
    0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	];

 xall = [];
 yall = [];
 colorAll = [];
 faceAlphaAll = [];
 counter = 1;

 for ii = 1:length(regions)

     thisRegion = combinedResults.(regions{ii}.label);
     color = cmap1(ii,:);
     faceAlpha = 0.8;

     x = ones(thisRegion.totalNumberOfCells,1)*(counter);
     y = thisRegion.faceExpPrefIndex';

     xall = cat(1,xall, x);
     yall = cat(1,yall, y);
     colorAll = cat(1,colorAll, color);
     faceAlphaAll = cat(1,faceAlphaAll, faceAlpha);
     counter = counter +2;

 end

plotFig = figure;
b = beeswarm(xall,yall,'Colormap', colorAll,'dot_size',1.2,'MarkerFaceAlpha',faceAlphaAll','overlay_style','ci','sort_style','fan');
xticks([1 3 5 7]);
xticklabels({'S1', 'M1','PMv', 'M3'})
xlabel('Region')
ylim([0 1.1]);
ylabel('Preference Index');
title('Face Expression Preference Indices per Regions')

 %% repeat for lateral v. medial areas
 % xall = [];
 %  yall = [];
 %  colorAll = [];
 %  counter = 1;
 % regions2plot{1}.label = 'lateralRegions';
 % regions2plot{2}.label = 'M3';
 %
 %
 % for ii = 1:length(regions2plot)
 %     for subj = 1:2
 %         if subj == 1
 %             thisRegion = barneyResults.(regions2plot{ii}.label);
 %             color = [0 0.4470 0.7410]; %barney blue
 %         elseif subj == 2
 %             thisRegion = thorResults.(regions2plot{ii}.label);
 %             color = [0.8500 0.3250 0.0980]; %thor orange
 %         end
 %         x = ones(thisRegion.totalNumberOfCells,1)*counter;
 %         y = thisRegion.faceExpPrefIndex';
 %
 %         xall = cat(1,xall, x);
 %         yall = cat(1,yall, y);
 %         colorAll = cat(1,colorAll, color);
 %         if subj == 1
 %             counter = counter +1;
 %         elseif subj == 2
 %             counter = counter +1.75;
 %         end
 %     end
 % end
 %
 % plotFig = figure;
 % b = beeswarm(xall,yall,'Colormap', colorAll,'dot_size',1.2);
 % xticks([1.5 4.25 ]);
 % xticklabels({'lateral cortex','medial cortex'})
 % xlabel('Region')
 % ylim([0 1.2]);
 % ylabel('faceExp Pref Index');
 % title('Face Expression Preference Indices per Regions')
 % legend('Barney', 'Thor');
 %
 %% one last time for M1 v. all others
 xall = [];
  yall = [];
  colorAll = [];
  counter = 1;
 regions2plot{1}.label = 'combinedM1';
 regions2plot{2}.label = 'nonM1';
 
 
for ii = 1:length(regions2plot)
    for subj = 1:2
        if subj == 1
            thisRegion = barneyResults.(regions2plot{ii}.label);
            color = [0 0.4470 0.7410]; %barney blue 
        elseif subj == 2
            thisRegion = thorResults.(regions2plot{ii}.label);
            color = [0.8500 0.3250 0.0980]; %thor orange 
        end
        x = ones(thisRegion.totalNumberOfCells,1)*counter;
        y = thisRegion.faceExpPrefIndex';

        xall = cat(1,xall, x);
        yall = cat(1,yall, y);
        colorAll = cat(1,colorAll, color);
        if subj == 1
            counter = counter +1;
        elseif subj == 2
            counter = counter +1.75;
        end
    end
end

plotFig = figure;
%b = beeswarm(xall,yall,'Colormap', colorAll,'sort_style','hex','dot_size',1.2);
b = beeswarm(xall,yall,'Colormap', colorAll,'dot_size',1.2,'MarkerFaceAlpha',faceAlphaAll','overlay_style','ci','sort_style','fan');

xticks([1.5 4.25 ]);
xticklabels({'all M1','PMv + M3'})
xlabel('Region')
ylim([0 1.1]);
ylabel('Preference Index');
title('Face Expression Preference Indices per Regions')
legend('Monkey1', '','','Monkey2');

