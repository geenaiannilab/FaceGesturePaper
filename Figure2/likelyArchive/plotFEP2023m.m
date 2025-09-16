%%% FEP = (nExp - [(sum (Ri) / rmax)]) ./ nExp -1 
%%% 
%%% Ri = spike rate for expression, i 
%%% Rmax = max Ri of any expression
%%% nExp = number of expressions compared 
%%%
%%% FEP = 0 indicates equal for all expressions
%%% FEP = 1 indicates discharge for 1 exp only 
%%% for building intuition -- for 3 expressions:
%%%  neuron firing 3hz, 1hz, 1hz --> 0.67
%%%  neuron firing 2hz, 1hz, 1hz --> 0.5
%%%  neuron firing 2hz, 2hz, 1hz --> 0.25


close all; clear all 
set(0,'defaultAxesFontSize', 36); % bc im blind 
set(0,'defaultAxesFontWeight', 'bold');

FEPthreshold = 0.5;
expColorMap = [1 0 0
    0 0 1
    0 1 0]; % r,b,g
scaleVal = 40; % for plotting
subjects = {'barney','thor'};
pathToMatFiles = '/Users/geena/Dropbox/PhD/SUAinfo/';

%% get barney's data 
barneyRawData1 = load(fullfile(pathToMatFiles, 'Barney_210704/Data4Analysis/210704_Barney_faceExpPref_V2.mat'));
barneyRawData2 = load(fullfile(pathToMatFiles, 'Barney_210706/Data4Analysis/210706_Barney_faceExpPref_V2.mat'));
barneyRawData3 = load(fullfile(pathToMatFiles, 'Barney_210805/Data4Analysis/210805_Barney_faceExpPref_V2.mat'));

barneyRawData.FEP = [barneyRawData1.faceExpPref, barneyRawData2.faceExpPref, barneyRawData3.faceExpPref]';
barneyRawData.meanFRs = cat(1, barneyRawData1.Ris, barneyRawData2.Ris, barneyRawData3.Ris);
barneyRawData.spikeLabels2plot = [barneyRawData1.spikeLabels2plot; barneyRawData2.spikeLabels2plot;  barneyRawData3.spikeLabels2plot]';
barneyRawData.cortRegion = [barneyRawData1.cortRegion barneyRawData2.cortRegion barneyRawData3.cortRegion]';

%% get thor's data 
thorRawData1 = load(fullfile(pathToMatFiles, 'Thor_171010/Data4Analysis/171010_Thor_faceExpPref_V2.mat'));
thorRawData2 = load(fullfile(pathToMatFiles, 'Thor_171027/Data4Analysis/171027_Thor_faceExpPref_V2.mat'));
thorRawData3 = load(fullfile(pathToMatFiles, 'Thor_171005/Data4Analysis/171005_Thor_faceExpPref_V2.mat'));
thorRawData4 = load(fullfile(pathToMatFiles, 'Thor_171128/Data4Analysis/171128_Thor_faceExpPref_V2.mat'));

thorRawData.FEP = [thorRawData1.faceExpPref, thorRawData2.faceExpPref, thorRawData3.faceExpPref, thorRawData4.faceExpPref]';
thorRawData.meanFRs = cat(1, thorRawData1.Ris, thorRawData2.Ris, thorRawData3.Ris, thorRawData4.Ris);
thorRawData.spikeLabels2plot = [thorRawData1.spikeLabels2plot; thorRawData2.spikeLabels2plot;  thorRawData3.spikeLabels2plot; thorRawData4.spikeLabels2plot]';
thorRawData.cortRegion = [thorRawData1.cortRegion thorRawData2.cortRegion thorRawData3.cortRegion thorRawData4.cortRegion]';

%% get the cortical regions per cell 
cortRegionOut = [];
FEPout = [];
meanFRsOut = [];

for ss = 1:length(subjects)
    if strcmpi(subjects{ss},'Barney')
        data = barneyRawData;
    elseif strcmpi(subjects{ss},'Thor')
        data = thorRawData;
    end
    
    % pool faceExpPrefs indices 
    FEP = data.FEP;
    cortRegion = data.cortRegion; 
    meanFRs = data.meanFRs;

    cortRegionOut = [cortRegionOut; cortRegion];
    FEPout = [FEPout ; FEP];
    meanFRsOut = [meanFRsOut; meanFRs];

    clear cortRegion; clear FEP; clear meanFRs;
end

cellsAboveFEPThresh = FEPout >= FEPthreshold;
meanFRs = meanFRsOut(cellsAboveFEPThresh,:);



xRegions = categorical(cortRegionOut(cellsAboveFEPThresh));
xRegions = reordercats(xRegions,{'S1' 'M1' 'PMv','M3'});
yFEPvalues = FEPout(cellsAboveFEPThresh);

[cExpWinnerFR, cExpWinner] = max(meanFRs,[],2);

%%
figure; 
%swarmchart(xRegions,yFEPvalues,cExpWinnerFR*scaleVal,cExpWinner,'filled','XJitter','density','XJitterWidth',1.1)
swarmchart(xRegions,yFEPvalues, 200, cExpWinner,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.5)%'XJitter','density','XJitterWidth',1.1);
colormap(expColorMap)
yticks(0.5:0.1:1)
ylim([.45 1.05])
ylabel('Preference Index')
xlabel('Region')
title('Preference Indices of Selective Cells')
axes = gca;
axes.FontWeight = 'bold';


