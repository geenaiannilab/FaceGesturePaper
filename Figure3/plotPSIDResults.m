%%% assumes you've run analyzePseudopopPSID2023.m
% Figures plotted 
%   Average CC performance per region (bhvPC 1-5)
% 	Distribution of optimal hyper parameters per region 
% 	Neural, full model self-prediction accuracy 
% 	Neural, first stage model self-prediction accuracy 
% 	(ratio of two) 
clear all; close all 

set(0,'defaultAxesFontSize', 24); % bc im blind 
set(0,'defaultAxesFontWeight', 'bold'); % bc im blind 

% define the dataset you are running
thisSubj = 'combined';
pseudoPop = 'LSOnly/balanced'; % vs 'unbalanced'
subj = 'combinedSubjs'; % vs 'separateSubjs' v 'combinedSubjs'
nIterations = 25; % pseduopopulations
popSize = '50'; % nCells

workDir = ['/Users/geena/Dropbox/PhD/SUAinfo/PSIDResults2023/Pseudopopulations/' pseudoPop '/' subj '/N_' num2str(popSize) 'cells/results'];

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

%  
results = load([workDir '/decodingResults.mat']);
regionList = fieldnames(results.resultsOut.combined)';

x = 2*[ones(1,nIterations) ones(1,nIterations)*2 ones(1,nIterations)*3 ones(1,nIterations)*4 ones(1,nIterations)*5] ;

allDecodingMean = [];
allDecodingErr = [];
alldistN1 = [];
alldistNX = [];
allIter2plot = [];

for rr = 1:length(regionList) 
    theseResults = results.resultsOut.combined.(regionList{rr});
    outcomeMetric = theseResults.outcomeMetric;

    %%%  kinematic decoding is nFolds x nBhvsPred x nIter
    [~, maxIdx2plot] = max(squeeze(mean(results.resultsOut.combined.(regionList{rr}).kinematicDecoding(:,1,:),1)));
    decodingAcc = round(squeeze(mean(theseResults.kinematicDecoding,1,'omitnan')),2); % averaged over nOuterFolds
    [decodingErr, decodingMean] = std(decodingAcc,[],2);
    
    allIter2plot(:,rr) = maxIdx2plot;
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

end

% plot decoding Means (nBhvsPCs x region) w/ 2std errorbars

figure;
b = bar(allDecodingMean);
for rr = 1:length(regionList)
    b(rr).FaceColor = regionCmap(rr,:);
end
hold on;

ngroups = size(allDecodingMean, 1); nbars = size(allDecodingMean, 2); groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    tmp = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(tmp, allDecodingMean(:,i), 2*allDecodingErr(:,i),'HandleVisibility','off');
    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 3;
end
hold off
legend(regionList); xlabel('Face Motion PC'); ylabel(outcomeMetric)
title(['Mean ' outcomeMetric ' : Decoding ' extractBefore(pseudoPop,'/') ])

%%
% N1 and Nx distributions
figure;
subplot 121
beeswarm(x',alldistN1','sort_style','rand','overlay_style','ci','colormap',regionCmap);
axes1 = gca; axes1.XTick = 2:2:2*length(regionList); axes1.XTickLabel = regionList; 
ylabel('N1: Optimal # Behaviorally Relevant States')

subplot 122
beeswarm(x',alldistNX','sort_style','rand','overlay_style','ci','colormap',regionCmap);
axes2 = gca; axes2.XTick = 2:2:2*length(regionList); axes2.XTickLabel = regionList; 
ylabel('NX: Optimal Total # States')

sgtitle(['Optimal N1 and NX : ' extractBefore(pseudoPop,'/')]);

% neural self prediction FULL model
figure;
subplot 121
x = categorical(regionList); x = reordercats(x, {'S1','M1','PMv','M3','All'});
b = bar(x, mean(allNeuralSelfFullMean,1),'facecolor','flat'); b.CData = regionCmap;
hold on; 

er = errorbar(x,mean(allNeuralSelfFullMean,1), mean(allNeuralSelfFullErr,1),'HandleVisibility','off');
er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 3;
hold off; 
ylabel(outcomeMetric);ylim([0 1]);
title(['Neural Self Predict, Full Model ' extractBefore(pseudoPop,'/')]);

% neural self prediction FIRST stage only 
subplot 122;
x = categorical(regionList); x = reordercats(x, {'S1','M1','PMv','M3','All'});
b = bar(x, mean(allNeuralSelfFirstMean,1),'facecolor','flat'); b.CData = regionCmap;
hold on; 

er = errorbar(x,mean(allNeuralSelfFirstMean,1), mean(allNeuralSelfFirstErr,1),'HandleVisibility','off');
er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 3;
hold off; 
ylabel(outcomeMetric); ylim([0 1]);
title(['Neural Self Predict, First Stage Only Model ' extractBefore(pseudoPop,'/')]);

%  disp(['Average Neural, Full Model ' outcomeMetric ' : ' num2str(meanSelfPred_2stageAcc) ' +/- ' num2str(stdSelfPred_2stageAcc)]);
%    disp(['Average Neural First-stage only ' outcomeMetric ' : ' num2str(meanSelfPred_1stageAcc) ' +/- ' num2str(stdSelfPred_1stageAcc)]);
%    disp(['Kinematic : AllOther Info Ratio = ' num2str((meanSelfPred_1stageAcc / meanSelfPred_2stageAcc))])

ratio = mean(allNeuralSelfFirstMean,1) ./ mean(allNeuralSelfFullMean,1);
ratioErr = [mean(allNeuralSelfFirstErr,1) ; mean(allNeuralSelfFullErr,1)];


%%%%%%%% 
%% plot iteration w/ best decoding of PC1

for rr = 1:length(regionList) 
    
    iter2plot = allIter2plot(:,rr);

    theseResults = results.resultsOut.combined.(regionList{rr});

    outerBhvTestC = cell2mat(theseResults.bhvTest{iter2plot});
    outerBhvPredC = cell2mat(theseResults.bhvPred{iter2plot});
    outerBhvTest = theseResults.bhvTest{iter2plot};
    
    heldOutBhvTest = reshape(outerBhvTestC,length(taxis),length(outerBhvTest),nz);
    heldOutBhvPred = reshape(outerBhvPredC,length(taxis),length(outerBhvTest),nz);
    
    figure; 
    ha = tight_subplot(2,2,0.125,[0.075 0.1],[0.1 0.05]);

    for bb = 1:nz-1 % each predicted bhvPC
        
        test1Y = mean(heldOutBhvTest(:,:,bb),2);
        test1Err = 2*(std(heldOutBhvTest(:,:,bb),[],2) ./ sqrt(length(outerBhvTest)));
    
        pred1Y = mean(heldOutBhvPred(:,:,bb),2);
        pred1Err = 2*(std(heldOutBhvPred(:,:,bb),[],2) ./ sqrt(length(outerBhvTest)));
        
        axes(ha(bb))
        shadedErrorBar(taxis,test1Y,test1Err,'lineprops',{'Color','k','linew',8});
        hold on;
        shadedErrorBar(taxis, pred1Y,pred1Err,'lineprops',{'Color',[0.3010 0.5 0.9],'linew',8});
        hold off

        title(['Face PC' num2str(bb)],'FontSize',26)
        set(ha(bb), 'YTickLabelMode','auto');
        xlim([-0.5 1.5])
        if bb == 1 
            xlabel('Time, s'); ylabel('a.u.'); 
            legend('true','pred')
        end

    end

    set(ha(1:nz-1),'XTickLabel',[-0.5 0 0.5 1 1.5]); 

    sgtitle([regionList{rr}],'FontSize',36,'FontWeight','bold')
    thisPlot = gcf; thisPlot.Position = [ 476    86   656   780];
    plotFig = gcf; 
    %exportgraphics(plotFig, ['~/Dropbox/subFigurePDFs/PSID_kinematics_' extractBefore(pseudoPop,'/') '_' regionList{rr} '.pdf'], 'ContentType', 'vector');
end

%%
for bb = 1:nz-1

    figure; 
    subplot 211
    for tt = 1:length(heldOutBhvTest)
        plot(taxis, heldOutBhvTest(:,tt,bb),'b'); hold on 

    end
    plot(taxis, mean(heldOutBhvTest(:,:,bb),2),'k','linew',4); hold off

    subplot 212 
    for tt = 1:length(heldOutBhvPred)
        plot(taxis, heldOutBhvPred(:,tt,bb),'r'); hold on 

    end
    plot(taxis, mean(heldOutBhvPred(:,:,bb),2),'k','linew',4); hold off

end

% neural trajectories 
colorTAxis = repmat(taxis,numel(outerLatentPredC(:,1)),1);
colorTAxis = colorTAxis(1,:); markerSize = 52;

grey = [179 179 179]/255 ;
red = [255 0 0]/255 ;
blue = [0 0 255]/255 ;
green = [0 204 0 ]/255;
black = [0 0 0 ]/255;

n = length(taxis) ; 
cmapBhv(:,:,1) = [linspace(grey(1),red(1),n)' linspace(grey(2),red(2),n)' linspace(grey(3),red(3),n)'];
cmapBhv(:,:,2) = [linspace(grey(1),blue(1),n)' linspace(grey(2),blue(2),n)' linspace(grey(3),blue(3),n)'];
cmapBhv(:,:,3) = [linspace(grey(1),green(1),n)' linspace(grey(2),green(2),n)' linspace(grey(3),green(3),n)'];
cmapBhv(:,:,4) = [linspace(grey(1),black(1),n)' linspace(grey(2),black(2),n)' linspace(grey(3),black(3),n)'];

%%% plot neural trajectories (per bhv)

for ii = 1:length(bhvs2test)
    figure;
    bhv = bhvs2test(ii);
    theseTrials = cell2mat(outerLatentPred(outerBhvTestLabels == bhv));
    theseTrials = reshape(theseTrials,length(taxis),length(outerLatentPred(outerBhvTestLabels == bhv)),optTotalStates);
    
    trialAvg = squeeze(mean(theseTrials,2));
    trialErr = 2*(squeeze((std(theseTrials,[],2))) ./ sqrt(length(outerLatentPred(outerBhvTestLabels == bhv))));
    
    subplot 211
    for jj = 1:optBhvRelevantStates
        plot(taxis, trialAvg(:,jj),'linew',3); hold on; 
    end
    hold off; legend; title('BhvRelev Dimensions')

    subplot 212
    for jj = optBhvRelevantStates+1:1:optTotalStates
        plot(taxis, trialAvg(:,jj),'linew',3); hold on; 
    end
    hold off; legend; title('BhvIrr Dimensions')
    sgtitle(['Bhv = ' num2str(bhvs2test(ii))])

    figure;
    scatter3(trialAvg(:,1),trialAvg(:,2),trialAvg(:,3),markerSize,cmapBhv(:,:,ii),'filled');
    title(['Bhv = ' num2str(bhvs2test(ii))])
end


%
figure; 
for ii = 1:length(bhvs2test)
    bhv = bhvs2test(ii);
    theseTrials = cell2mat(outerLatentPred(outerBhvTestLabels == bhv));
    theseTrials = reshape(theseTrials,length(taxis),length(outerLatentPred(outerBhvTestLabels == bhv)),optTotalStates);
    
    trialAvg = squeeze(mean(theseTrials,2));
    trialErr = 2*(squeeze((std(theseTrials,[],2))) ./ sqrt(length(outerLatentPred(outerBhvTestLabels == bhv))));
    
   scatter3(trialAvg(:,1),trialAvg(:,2),trialAvg(:,3),markerSize,cmapBhv(:,:,ii),'filled'); hold on; 
 
end
%%

