clear all; 
set(0,'defaultAxesFontSize',20)

%% intial set-up 
subject = 'combined';
nIterations = 50;
windowLength = 0.4;
windowShiftSize = 0.05;

workdir = ['/Users/geena/Dropbox/PhD/SUAinfo/Pseudopopulations/balanced/N_50cells/SVM_' subject '/'];

for iter = 1:nIterations 
    results(iter) = load([workdir '/' subject '_N' num2str(iter) '_SVMresults_' num2str(windowLength) '_' num2str(windowShiftSize) '.mat']);
    perm(iter) = load([workdir '/' subject '_N' num2str(iter) '_permSVMresults_' num2str(windowLength) '_' num2str(windowShiftSize) '.mat']);
end

%% define time axis 
firstWindowCenter = results(1).desiredStart + results(1).windowLength/2;
lastWindowCenter = results(1).desiredEnd - results(1).windowLength/2;
windowsCenter = firstWindowCenter:results(1).windowShiftSize:lastWindowCenter;

arrayList = fieldnames(results(1).validationScoreStruct);
validationScores = [results(:).validationScoreStruct];
SEMScores = [results(:).SEM];
%STDScores = [results(:).STD];
permScores = [perm(:).permValidationScoreStruct];

%% get averages across N pseudopopulations
for aa = 1:length(arrayList)

    allIterResults = [];
    allIterSEMs = [];

    for iter = 1:nIterations

        thisIterResults = validationScores(iter).(arrayList{aa});
        allIterResults = [thisIterResults; allIterResults];

        thisIterSEM = SEMScores(iter).(arrayList{aa});
        allIterSEMs = [thisIterSEM; allIterSEMs];
        
    end

    meanIterValidationScores(aa,:) = mean(allIterResults,1);
    stdErrBars(aa,:) = std(allIterResults,[],1); 
    semErrBars(aa,:) = mean(allIterSEMs,1);
    semErrorBars_2(aa,:) = std(allIterResults,1) ./ sqrt(nIterations) ;

end

%% permutations 98th% as cutoff for significance of decoder 
prcThres = 98;
P = zeros(length(arrayList),length(windowsCenter),nIterations); 

for aa = 1:length(arrayList)

    for iter = 1:nIterations

        for tt = 1:length(windowsCenter)

            P(aa,tt,iter) = prctile(permScores(iter).(arrayList{aa})(tt,:),prcThres);

        end

    end

end

%% plot the results

figure;

colorArray = [0.4940 0.1840 0.5560;0.6350 0.0780 0.1840;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];

for aa = 1:length(arrayList)
  
    chance = mean(P(aa,:,:),3)*100;
    chanceSEM = (std(P(aa,:,:),[],3)./ sqrt(nIterations) ) *100;
    
    errorbar(windowsCenter, meanIterValidationScores(aa,:)*100, stdErrBars(aa,:)*100,...
        'Color', colorArray(aa,:), 'LineWidth',8); hold on
    errorbar(windowsCenter,chance, chanceSEM,'Color', [0.5 0.5 0.5],'LineWidth',3); hold on
  
end
hold off
xlim([- 0.8 0.8]);
axes = gca;
axes.FontWeight = 'bold'
axes.FontSize = 34;
axes.XTick = windowsCenter(1:4:end);
xlabel('Time (s)')
ylabel('Mean Accuracy (%)')
ylim([50 103])
title(['Decoding Performance By Region'],'FontSize',45);

%%
figure;

colorArray = [0.4940 0.1840 0.5560;0.6350 0.0780 0.1840;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];

for aa = 1:length(arrayList)
  
    chance = mean(P(aa,:,:),3)*100;
    chanceSEM = (std(P(aa,:,:),[],3)./ sqrt(nIterations) ) *100;
    
    errorbar(windowsCenter, meanIterValidationScores(aa,:)*100, semErrBars(aa,:)*100,...
        'Color', colorArray(aa,:), 'LineWidth',8); hold on
    errorbar(windowsCenter,chance, chanceSEM,'Color', [0.5 0.5 0.5],'LineWidth',3); hold on
  
end
hold off
xlim([- 0.8 0.8]);
axes = gca;
axes.FontWeight = 'bold'
axes.FontSize = 34;
axes.XTick = windowsCenter(1:4:end);
xlabel('Time (s)')
ylabel('Mean Accuracy (%)')
ylim([50 103])
title(['Decoding Performance By Region #2'],'FontSize',45);

%%%
figure;

colorArray = [0.4940 0.1840 0.5560;0.6350 0.0780 0.1840;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];

for aa = 1:length(arrayList)
  
    chance = mean(P(aa,:,:),3)*100;
    chanceSEM = (std(P(aa,:,:),[],3)./ sqrt(nIterations) ) *100;
    
    errorbar(windowsCenter, meanIterValidationScores(aa,:)*100, semErrorBars_2(aa,:)*100,...
        'Color', colorArray(aa,:), 'LineWidth',8); hold on
    errorbar(windowsCenter,chance, chanceSEM,'Color', [0.5 0.5 0.5],'LineWidth',3); hold on
  
end
hold off
xlim([- 0.8 0.8]);
axes = gca;
axes.FontWeight = 'bold'
axes.FontSize = 34;
axes.XTick = windowsCenter(1:4:end);
xlabel('Time (s)')
ylabel('Mean Accuracy (%)')
ylim([50 103])
title(['Decoding Performance By Region #3'],'FontSize',45);


figure;
for aa = 1:length(arrayList)
    subplot(2,2,aa)
    errorbar(windowsCenter, meanIterValidationScores(aa,:)*100, semIterValidationScores(aa,:)*100,'LineWidth',8); hold on
    
    %for iter = 1:nIterations
    chance = mean(P(aa,:,:),3)*100;
    chanceSEM = (std(P(aa,:,:),[],3)./ sqrt(nIterations) ) *100;
    errorbar(windowsCenter,chance, chanceSEM,'LineWidth',3); hold on
    %end
    %hold off; 
    xlim([-0.8 0.8])
    xlabel('Time (s)')
    ylabel('Mean Decoding Accuracy (%)')
    ylim([50 105])
    title([arrayList{aa}])
    sgtitle([subject '- Pseudopopulations, N= ' num2str(nIterations)],'FontSize',20);

end
%%

% plot per iteration result 
figure;

for iter = 1:50
    for aa = 1:length(arrayList)
        subplot(2,2,aa)
        errorbar(windowsCenter, results(iter).validationScoreStruct.(arrayList{aa})*100, results(iter).SEM.(arrayList{aa})*100,'LineWidth',8); hold on

        plot(windowsCenter,P(aa,:,iter)*100,'LineWidth',3); hold on
        hold off;
        xlim([-0.8 0.8])
        xlabel('Time (s)')
        ylabel('Mean Decoding Accuracy (%)')
        ylim([40 110])
        title([arrayList{aa}])
        sgtitle([subject '- Pseudopopulations, N= ' num2str(nIterations)],'FontSize',20);
 
    end
    pause; clf
end


%
% figure; 
% 
% 
% for iter = 1:nIterations
% 
%     subplot 131; imagesc(results(iter).betaScoreStruct.S1(:,:,1)'); colormap(jet); caxis([-0.4 0.4])
%     title('Hyperplane, Thr vs LS');
%     
%     subplot 132; imagesc(results(iter).betaScoreStruct.S1(:,:,2)'); colormap(jet); caxis([-0.4 0.4])
%     title('Hyperplane, Thr vs Chew');
%     
%     subplot 133; imagesc(results(iter).betaScoreStruct.S1(:,:,3)'); colormap(jet);caxis([-0.4 0.4])
%     title('Hyperplane, LS vs Chew');
%     pause; clf
% end
% 


