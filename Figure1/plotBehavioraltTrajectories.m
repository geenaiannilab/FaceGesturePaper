
%%%%%%%%
%%%%%%%% PLOTTING Figure 1 
%%%%%%%%%
clear all; close all;
subject = 'Thor'; 
date = '171010';
workDir = ['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/DLC'];
load([workDir '/markerTrajectories.mat']);
scalar = 3000; 
tmin = taxis(1); tmax = taxis(end);
clear markerData2plot; 

% eigenspectrum
fig1 = figure; plot(cumsum(latent) ./sum(latent),'o','linew',8);
xlabel('PC'); ylabel('CumVar, %')
xlim([0 length(latent)])
title(['Marker Positions; Eigenspectrum; tmin = ' num2str(tmin) ' tmax = ' num2str(tmax)]);

% raster of PC weights
%fig2 = figure('Position', get(0, 'Screensize'));
%ax = imagesc(coeff(:,1:10)); colormap(bluewhitered);
%yticks = 1:2:size(allBhv,2);
%yticklabels = markerLabels;
%set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
%xlabel('PC'); title(['Weights, Markers tmin = ' num2str(tmin) ' tmax = ' num2str(tmax)])
%saveas(fig, ['~/Desktop/bhvWeights_tmin' num2str(tmin)  '.png']);

%% 2D depiction of PC weights
for PC = 1:5
    fig3 = figure('Position', [ 668    86   845   780]);
    x=[]; y = []; radius = [];
    for i = 1:length(markerLabels)
        xpos = (i*2) - 1;
        ypos = (i*2);
        x(i) = meanPos(xpos);
        y(i) = -meanPos(ypos);
        %radius(i) = mean([abs(PCweights(xpos,PC)) abs(PCweights(ypos,PC))]);
        radius(i) = PCweights(ypos,PC);
    end
    s = scatter(x,y,abs(radius)*scalar,radius,'filled'); 
    hold off; colormap(bluewhitered); axes = gca; axes.XTick = []; axes.XTickLabels = '';axes.YTick = []; axes.YTickLabels = '';
    title(['PC' num2str(PC) ' : ' num2str(explained(PC)) '% variance'])

end

%% trial-average & single-trial trajectories 
fig = figure('Position', get(0, 'Screensize')); 
set(gcf,'color','w');
colors = 'rbg';
markerType = 'osd*.'; markerSize = 300;

for bhv = 1:length(bhvTraj)

    data = (bhvTraj(bhv).data);

    for tt = 1:size(data,2)
        p = plot3(data(:,tt,1), data(:,tt,2),data(:,tt,3),...
            'Color',colors(bhv),'linew',0.5); hold on;
    end

end

for bhv = 1:length(bhvTraj)
    scatter3(bhvTraj(bhv).avg(:,1), bhvTraj(bhv).avg(:,2),bhvTraj(bhv).avg(:,3),...
        markerSize, colormap(customColormap(:,:,bhv)),'filled',markerType(bhv));  hold on     
end

hold off;
xlabel('Dim 1; au'); 
ylabel('Dim 2; au');
zlabel('Dim 3; au');
title(['Bhv Trajectories; tmin = ' num2str(tmin) ' tmax = ' num2str(tmax)])
grid on

view(3)
set(0,'defaultAxesFontSize',20);

%% plot the PCs with mean/SEM
fig = figure('Position', get(0, 'Screensize')); 
for bhv = 1:length(bhvs2plot)
    for PC = 1:5
        subplot(1,5,PC)
        colors = 'rbmgc';

        shadedErrorBar(taxis,bhvTraj(bhv).avg(:,PC),bhvTraj(bhv).sem(:,PC),'lineprops',{'Color',colors(bhvs2plot(bhv)),'linew',8}); hold on;

        %  plot(downsampledTimeAxis,bhvTraj(1).avg(:,PC),'r','linew',6); hold on; ...
      %      plot(downsampledTimeAxis,bhvTraj(2).avg(:,PC),'b','linew',6)
      %      plot(downsampledTimeAxis,bhvTraj(3).avg(:,PC),'g','linew',6); hold off
        title(['PC=' num2str(PC)],'FontSize',20)
    end
end

legend('Thr','LS', 'Chew')
sgtitle(['Top 10 FaceMotion PCs; tmin = ' num2str(tmin) ' tmax = ' num2str(tmax)],'FontSize',24)

%% plot per trial PCs 
for bhv = 1:length(bhvs2plot)
    fig = figure('Position', get(0, 'Screensize')); 

    trials = size(bhvTraj(bhv).data,2);
    
    for PC = 1:3
    subplot(2,3,PC)

        for tt = 1:trials
            markerPos = squeeze(bhvTraj(bhv).data(:,tt,PC));
            plot(taxis, markerPos,'Color',colors(bhvs2plot(bhv))); hold on
        end 
        hold off ; ax = gca; ax.YAxis.Visible = 'off';
        title(['MarkerPosition, PC=' num2str(PC)],'FontSize',20)

        subplot(2,3,PC+3)
        for tt = 1:trials
            markerPos = squeeze(bhvTraj(bhv).data(:,tt,PC));
            markerVel = diff(markerPos) ./ diff(taxis);
            plot(taxis(1:end-1), markerVel,'Color',colors(bhvs2plot(bhv))); hold on ;
        end
        hold off; ax = gca; ax.YAxis.Visible = 'off';
        title(['MarkerVelocity, PC=' num2str(PC)],'FontSize',20)

    end
    sgtitle(['Bhv = ' gestureNames{bhv}],'FontSize',30);
end

%%
figure;
colors = 'rbg';
sampleFrac = 0.01;
nTotalTrials = sum([size(markerData2plot(1).data,2),size(markerData2plot(2).data,2),size(markerData2plot(3).data,2)]);
randSample = ceil(sampleFrac*nTotalTrials);

for bhv = 1:length(bhvs2plot)
    
    thisData = markerData2plot(bhv).data;
    
    for tt = 1:size(thisData,2)
        for i = 1:length(markerLabels)
            xpos = (i*2) - 1;
            ypos = (i*2);
            x = markerData2plot(bhv).data(1,tt,xpos);
            y = -markerData2plot(bhv).data(1,tt,ypos);

            s = scatter(x,y,150,colors(bhv),'filled'); 
            s.MarkerFaceAlpha = 0.1; hold on; 

        end
    end
end
title('Starting Position for All Trials');



% 
% figure;
% colors = 'rbg';
% 
% for bhv = 1:length(bhvs2plot)
%     
%     for i = 1:length(markerLabels)
%         xpos = (i*2) - 1;
%         ypos = (i*2);
%         x = markerData2plot(bhv).avg(:,xpos);
%         y = -markerData2plot(bhv).avg(:,ypos);
%         xErr = markerData2plot(bhv).sem(:,xpos)';
%         yErr = markerData2plot(bhv).sem(:,ypos)';
%         err = [xErr; yErr];
%         shadedErrorBar(x,y,err,'lineprops',{colors(bhv),'linew',5}); hold on
% 
%     end
% end
% title('Mean and SEM positions, all Trials')
% 
% 
% 
