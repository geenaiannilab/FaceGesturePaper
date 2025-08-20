
%%% reminder to run for files V2 (+/- 1 s for FEPI) vs no V2 (+/- 0.5 s for
%%% FEPI) -- same results 

clear all; 
set(0,'defaultAxesFontSize', 16); % bc im blind 
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

for ss = 1%:length(subjects)
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

% run 1 way anova for effect of cortical region on FEP 
% and run post-hoc tests 
[P_oneway,ANOVATAB_oneway,STATS_oneway] = anova1(FEPout, cortRegionOut,'off');
[c_oneway,m_oneway] = multcompare(STATS_oneway);

% run 2 way anova for effect of cortical region, expresison type, on mean
% FR & run post-hoc tests

inputFR = [meanFRsOut(:,1) ; meanFRsOut(:,2); meanFRsOut(:,3)];
T    = cell(1, size(meanFRsOut,1)); T(:) = {'Threat'}; L = cell(1, size(meanFRsOut,1)); L(:) = {'Lipsmack'}; C = cell(1, size(meanFRsOut,1)); C(:) = {'Chew'};
inputExp = [T'; L'; C'];
inputCortRegion = [cortRegionOut ; cortRegionOut; cortRegionOut];
[P_twoway,TABLE_twoway,STATS_twoway] = anovan(inputFR,{inputExp inputCortRegion},'model','interaction','varnames',{'inputExp','inputCortRegion'});
[c_twoway,m_twoway, ~, names] = multcompare(STATS_twoway,'Dimension', [1 2]);

% plot it 
factors = STATS_twoway.varnames;
groupNames = STATS_twoway.grpnames;
cmap = [0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	];
start = 1; stop = 3;
for rr = 1:length(groupNames{2,1}) % per region
    
    means = m_twoway(start:stop,1);
    errors = m_twoway(start:stop,2);
    
    b(rr) = bar(start:1:stop, means); 
    b(rr).FaceColor = cmap(rr,:);
    hold on

    er = errorbar(start:stop, means,errors,'HandleVisibility','off');
    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 4;

    start = start +3; stop = stop + 3;
end
hold off

tickLabels = cell(1, length(m_twoway)); tickLabels(:) = {'Threat','Lipsmack','Chew','Threat','Lipsmack','Chew','Threat','Lipsmack','Chew','Threat','Lipsmack','Chew'};
set(gca,'XTick',1:length(m_twoway), 'XTickLabel',tickLabels);
legend({'S1','M1','PMv','M3'});

title('MeanFR by Region & GestureType')
ylabel('mean FR');

%% plot again just colored differently
cmap2 = [1 0 0
    0 0 1
    0 1 0];
start = 1; skip = 3;
meansAll = [];
errorsAll = [];
for exp = 1:length(groupNames{1,1}) % per region
    max = length(m_twoway);

    means = m_twoway(start:skip:max,1);
    errors = m_twoway(start:skip:max,2);

    meansAll = [meansAll means];
    errorsAll = [errorsAll errors];
    start = start +1;

end
%%
figure; 
b = bar(meansAll); hold on;

ngroups = size(meansAll, 1);
nbars = size(meansAll, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, meansAll(:,i), errorsAll(:,i),'HandleVisibility','off');
    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 4;
end
hold off

b(1).FaceColor = cmap2(1,:);b(2).FaceColor = cmap2(2,:);b(3).FaceColor = cmap2(3,:);

tickLabels = cell(1, 4); tickLabels(:) = {'S1','M1','PMv','M3'};
set(gca,'XTick',1:length(m_twoway), 'XTickLabel',tickLabels);
legend({'Threat','Lipsmack','Chew'});

title('MeanFR by Region & GestureType')
ylabel('mean FR');

