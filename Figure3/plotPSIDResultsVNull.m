%%%% written 250201 
%%%% plots the output of PSID decoding (all regions; Figure 3A from
%%%% manuscript) 
%%%% Versus decoding performance of a shuffled, null dataset 

clear all; close all; 
set(0,'defaultAxesFontSize', 16); % bc im blind 


% define the dataset you are running
thisSubj = 'combined';
pseudoPop = 'LSonly/balanced'; % vs 'unbalanced'
subj = 'combinedSubjs'; % vs 'separateSubjs' v 'combinedSubjs'
nullModelFlag = true;  nullModelType = 'shuffleBhvTimepoints';
nNullShuffles = 50;

nPseudopops = 25; % pseduopopulations
popSize = '50'; % nCells
outcomeMetric = 'CC'; % vs R2 vs RSME etc
saveFlag = true; 

dataDir = (['/Users/geena/Dropbox/PhD/SUAinfo/PSIDResults2023/Pseudopopulations/' pseudoPop '/' subj '/N_' popSize 'cells/']);
resultsDir = ['/Users/geena/Dropbox/PhD/SUAinfo/PSIDResults2023/Pseudopopulations/' pseudoPop '/' subj '/N_' popSize 'cells/results'];

% load one pseudopop true and Null data
decodingResults = load([resultsDir '/decodingResults.mat' ]); %true results 
decodingShuffleResults =  load([resultsDir '/decodingNullResults_pp_N23.mat' ]); % dummy null data with first 2 shuffles; add to shuffles N3-N25

% define the dataset you are running
thisSubj = 'combined';
pseudoPop = 'LSOnly/balanced'; % vs 'unbalanced'
subj = 'combinedSubjs'; % vs 'separateSubjs' v 'combinedSubjs'
nIterations = 25; % pseduopopulations
popSize = '50'; % nCells
regionList = fieldnames(decodingResults.resultsOut.combined)';

% some info on the data 
outerFolds = 4; 
statesBhvRelv = [2 3 4 8 10 12 20]; %  N1 gridsearch  
additionalStates = [2 4 8 10]; % ** this is NOT Nx (Nx = statesBhvRelv + additional) **
horizonValues = [5 10 15  ];  % i gridsearch 
nz = 5;                       % bhvDims predicted      
taxis = -0.5:0.05:1.5; % didnt save in results file, ugh

% region colormap
regionCmap = [0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250
    0.5 0.5 0.5];

x = 2*[ones(1,nIterations) ones(1,nIterations)*2 ones(1,nIterations)*3 ones(1,nIterations)*4 ones(1,nIterations)*5] ;

allDecodingMean = [];
allDecodingErr = [];
alldistN1 = [];
alldistNX = [];
allIter2plot = [];

for rr = 1:length(regionList) 
    theseResults = decodingResults.resultsOut.combined.(regionList{rr});
    shuffledResults = decodingShuffleResults.shuffledResultsAll.combined.(regionList{rr});

    outcomeMetric = theseResults.outcomeMetric;

    %%%  kinematic decoding is nFolds x nBhvsPred x nIter
    decodingAcc = round(squeeze(mean(theseResults.kinematicDecoding,1,'omitnan')),2); % averaged over nOuterFolds
    [decodingErr, decodingMean] = std(decodingAcc,[],2);
    
    % nBhvPreds x nRegions
    allDecodingMean(:,rr) = round(decodingMean,2);
    allDecodingErr(:,rr) = round(decodingErr,2) ./ sqrt(nIterations);

    %%% optimal state #s
    alldistN1 = [alldistN1 theseResults.N1];
    alldistNX = [alldistNX theseResults.NX ];

    %%% Neural, full model self-prediction accuracy
    neuralSelfFullAcc = round(squeeze(mean(theseResults.neuralSelfPred2stage,1,'omitnan')),2);
    [neuralSelfFullErr, neuralSelfFullMean] = std(neuralSelfFullAcc,[],2);
    allNeuralSelfFullMean(:,rr) = round(neuralSelfFullMean,2);
    allNeuralSelfFullErr(:,rr) = round(neuralSelfFullErr,2);

    neuralSelfFirstAcc = round(squeeze(mean(theseResults.neuralSelfPred1stage,1,'omitnan')),2);
    [neuralSelfFirstErr, neuralSelfFirstMean] = std(neuralSelfFirstAcc,[],2);
    allNeuralSelfFirstMean(:,rr) = round(neuralSelfFirstMean,2);
    allNeuralSelfFirstErr(:,rr) = round(neuralSelfFirstErr,2);

    % shuffled kinematic decoding is nBhvsPred x nIter
    % allShuffledMean is nBhvsPred x nIter x nRegions
    %allShuffledMean(:,:,rr) = round(shuffledResults.kinematicDecoding,2);
    %allShuffledStd(:,:,rr) = round(shuffledResults.kinematicDecodingStd,2);
    allShuffledMean(:,rr) = mean(shuffledResults.kinematicDecoding, 2);
    allShuffledStd(:,rr) = mean(shuffledResults.kinematicDecodingStd, 2);
    allShuffledErr(:,rr) = allShuffledStd(:,rr); % I'm not sure what error you want to plot!



    % these are all nCells x nIter x nBhvPred
    allNeuralSelfFullMeanShuffled(:,:,rr) = round(shuffledResults.neuralSelfPred2stage,2);
    allNeuralSelfFullMeanShuffledErr(:,:,rr) = round(shuffledResults.neuralSelfPred2stageStd,2);

    allNeuralSelf1stageMeanShuffled(:,:,rr) = round(shuffledResults.neuralSelfPred1stage,2);
    allNeuralSelf1stageMeanShuffledErr(:,:,rr) = round(shuffledResults.neuralSelfPred1stageStd,2);

end

%%
figure;
b = bar(allDecodingMean);
for rr = 1:length(regionList)
    b(rr).FaceColor = regionCmap(rr,:);
    b(rr).FaceAlpha = 0.5;
end
hold on;

ngroups = size(allDecodingMean, 1); nbars = size(allDecodingMean, 2); groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    tmp = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(tmp, allDecodingMean(:,i), 2*allDecodingErr(:,i),'HandleVisibility','off');
    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 1.5;
end

% plot the shuffled decoding accuracies 
b2 = bar(squeeze(allShuffledMean(:,:))); hold off
for rr = 1:length(regionList)
    b2(rr).FaceColor = regionCmap(rr,:);
    b2(rr).BarWidth = 0.4;
    b2(rr).FaceAlpha = 0.75;
end
hold on;

for i = 1:nbars
    tmp = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(tmp, allShuffledMean(:,i), 2*allShuffledErr(:,i),'HandleVisibility','off');
    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 1.5;
end

%legend(regionList); xlabel('Face Motion PC'); ylabel(outcomeMetric)
%title(['Mean ' outcomeMetric ' : Decoding ' extractBefore(pseudoPop,'/') ])
%

ylim([-0.1 0.8])
yticks([])
set(gcf, 'PaperSize', [5 5]); % Adjust these values based on your figure dimensions
set(gcf, 'PaperPosition', [0 0 5 5]);
print(gcf, '~/Desktop/myFigure.pdf', '-dpdf');