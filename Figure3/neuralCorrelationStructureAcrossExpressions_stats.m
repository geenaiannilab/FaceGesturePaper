% written GI 
%%%     
%%% PUPRPOSE: save neural signal correlation between all neuron pairs over time, per behavior (trial-avgeraged) 
%%%                  over ALL DAYS (to eventually get the mean) 

close all; clear all
set(0,'defaultAxesFontSize', 36); % bc im blind 
set(0,'defaultAxesFontWeight', 'bold'); 

dates = {'210704','210706','210805'};
subject ='Barney';
regions = {'S1','M1','PMv','M3','All'};
bhvs2plot = [1 2 4];
bhvStrings = {'Threat','Lipsmack','Chew'};

for dd = 1:length(dates)
    
    workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' dates{dd} '/CorrMatricesV2']);
    disp(['Date is ' dates{dd}]);

    for rr = 1:length(regions)
        data = load([workdir '/' regions{rr} '.mat']);
        
        disp(['Region is ' regions{rr}])

        taxis = data.taxis;

        % plot pairwise correlations
        tmp1 = fitlm(data.prCorr{:,1}', data.prCorr{:,2}');
        r2_1v2(dd,rr) = tmp1.Rsquared.Ordinary;
        p_1v2(dd,rr) = tmp1.ModelFitVsNullModel.Pvalue;

        tmp2 = fitlm(data.prCorr{:,1}', data.prCorr{:,3}');
        r2_1v3(dd,rr) = tmp2.Rsquared.Ordinary;
        p_1v3(dd,rr) = tmp2.ModelFitVsNullModel.Pvalue;

        tmp3 = fitlm(data.prCorr{:,2}', data.prCorr{:,3}');
        r2_2v3(dd,rr) = tmp3.Rsquared.Ordinary;
        p_2v3(dd,rr) = tmp3.ModelFitVsNullModel.Pvalue;

        clear tmp1 tmp2 tmp3;
    end
end

save(['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_R2ValsAcrossDays.mat'], 'dates','subject', 'regions', 'r2_1v3', 'r2_1v2', 'r2_2v3','p_1v2','p_1v3','p_2v3');

cmap1 = [         
    0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	
    0.75 0.75 0.75     ];

%%
%%% after you run the above for both subjects 
thorData = load('/Users/geena/Dropbox/PhD/SUAinfo/Thor_R2ValsAcrossDays.mat');
barneyData = load('/Users/geena/Dropbox/PhD/SUAinfo/Barney_R2ValsAcrossDays.mat');

thorMeans_1v2 = mean(thorData.r2_1v2,1); 
thorMeans_1v3 = mean(thorData.r2_1v3,1);
thorMeans_2v3 = mean(thorData.r2_2v3,1);

barneyMeans_1v2 = mean(barneyData.r2_1v2,1);
barneyMeans_1v3 = mean(barneyData.r2_1v3,1);
barneyMeans_2v3 = mean(barneyData.r2_2v3,1);

combinedData_1v2 = [thorData.r2_1v2; barneyData.r2_1v2];
combinedData_1v3 = [thorData.r2_1v3; barneyData.r2_1v3];
combinedData_2v3 = [thorData.r2_2v3; barneyData.r2_2v3];

combinedPval_1v2 = [thorData.p_1v2; barneyData.p_1v2];
combinedPval_1v3 = [thorData.p_1v3; barneyData.p_1v3];
combinedPval_2v3 = [thorData.p_2v3; barneyData.p_2v3];


combinedMeans_1v2 = mean([thorData.r2_1v2 ; barneyData.r2_1v2],1);
combinedMeans_1v3 = mean([thorData.r2_1v3 ; barneyData.r2_1v3],1);
combinedMeans_2v3 = mean([thorData.r2_2v3 ; barneyData.r2_2v3],1);

combinedStd_1v2 = std([thorData.r2_1v2 ; barneyData.r2_1v2],[],1);
combinedStd_1v3 = std([thorData.r2_1v3 ; barneyData.r2_1v3],[],1);
combinedStd_2v3 = std([thorData.r2_2v3 ; barneyData.r2_2v3],[],1);

combinedMeansAllBhv = mean([thorData.r2_1v2 ; thorData.r2_1v3; thorData.r2_2v3; barneyData.r2_1v2 ;barneyData.r2_1v3; barneyData.r2_2v3 ],1);
combinedStdAllBhv = std([thorData.r2_1v2 ; thorData.r2_1v3; thorData.r2_2v3; barneyData.r2_1v2 ;barneyData.r2_1v3; barneyData.r2_2v3 ],[],1);

combinedMeansAllBhvAllRegions = mean([thorData.r2_1v2 ; thorData.r2_1v3; thorData.r2_2v3; barneyData.r2_1v2 ;barneyData.r2_1v3; barneyData.r2_2v3 ],'all');
combinedStdAllBhvAllRegions = std([thorData.r2_1v2 ; thorData.r2_1v3; thorData.r2_2v3; barneyData.r2_1v2 ;barneyData.r2_1v3; barneyData.r2_2v3 ],[],'all');


% Plot R2 for threat v lipsmack (bhv1v2)
figure; 
counter = 1;

for dd = 1:size(combinedData_1v2,1)
    for rr = 1:length(regions)
        
        
        jitter = rand(1,1)*0.5;
        if rem(counter, 2) == 0
            x = rr*2 + jitter;
        else
             x = rr*2 - jitter;
        end
        s = swarmchart(x, combinedData_1v2(dd,rr),250,cmap1(rr,:),'filled'); hold on;
        alpha = (1-combinedPval_1v2(dd,rr))+0.1;
        if alpha > 1
            s.MarkerFaceAlpha = 1;
        else
            s.MarkerFaceAlpha = alpha;
        end
        
        counter = counter+1;

    end
end

% Plot R2 for threat V chew (bhv1v3)
counter = 1;

for dd = 1:size(combinedData_1v3,1)
    for rr = 1:length(regions)
        
        
        jitter = rand(1,1)*0.5;
        if rem(counter, 2) == 0
            x = rr*2 + jitter;
        else
             x = rr*2 - jitter;
        end
        s = swarmchart(x, combinedData_1v3(dd,rr),250,cmap1(rr,:),'filled'); hold on;
        alpha = (1-combinedPval_1v3(dd,rr))+0.1;
        if alpha > 1
            s.MarkerFaceAlpha = 1;
        else
            s.MarkerFaceAlpha = alpha;
        end
        
        counter = counter+1;

    end
end

% Plot R2 for lipsmack V chew (bhv2v3)

counter = 1;

for dd = 1:size(combinedData_2v3,1)
    for rr = 1:length(regions)
        
        
        jitter = rand(1,1)*0.5;
        if rem(counter, 2) == 0
            x = rr*2 + jitter;
        else
             x = rr*2 - jitter;
        end
        s = swarmchart(x, combinedData_2v3(dd,rr),250,cmap1(rr,:),'filled'); hold on;
        alpha = (1-combinedPval_2v3(dd,rr))+0.1;
        if alpha > 1
            s.MarkerFaceAlpha = 1;
        else
            s.MarkerFaceAlpha = alpha;
        end
        
        counter = counter+1;

    end
end

for rr = 1:length(regions)
    scatter(rr*2, combinedMeansAllBhv(rr),500,cmap1(rr,:),'+','linew',4);
    hold on 
end
hold off
ylabel('R-Squared');
title('Pairwise Neuron Correlations');

hold off; 
axes = gca;
axes.XTick = 2:2:10;
axes.XTickLabel = regions;

%%
