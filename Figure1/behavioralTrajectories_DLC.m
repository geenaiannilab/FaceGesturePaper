%% updated July 2023 to incorporate Thor's data, and do trial-chunking BEFORE decomposition 
%% tested and runs with old barney data too 

close all; clear all;
set(0,'defaultAxesFontSize',24);
set(0,'defaultAxesFontWeight','bold');

%% data params 
subject = 'Barney';
date = '210704';
runTSNE = false;
subsessions2plot = [2:5 7:10]; %
bhvs2plot = [1 2 4];
markers2use =  'all'; excludeTongue = 1; plotDLCFlag = false; 
newBinSize = 0.01;
tmin = 1;
tmax = 1;
centerNormalizeFlag = 0;
centerFlag = 1;
smoothSize = 10; 

minRestFlag = 1;
minRest = abs(tmin); % in sec; minimal rest prior to move onset (trials to include)

%% set up files to load 

workDir = ['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis'];
infoDir = ([workDir 'TDTdata/Info']);

%% concatenate subsession files 
% get relevant files 
flist = natsortfiles(dir([workDir '/bhvspikes_sub*.mat']));
flist = fullfile({flist.folder}, {flist.name}); 

% determine which type of subsession to plot
if strcmpi(subsessions2plot,'face expression')
    subs2use = flist(contains(flist, 'face expression'));
elseif strcmpi(subsessions2plot,'visual expression')
    subs2use = flist(contains(flist, 'visual expression'));
elseif strcmpi(subsessions2plot,'chew expression')
    subs2use = flist(contains(flist, 'chew expression'));
elseif strcmpi(subsessions2plot,'chew')
    subs2use = flist(contains(flist, 'chew'));
elseif strcmpi(subsessions2plot,'rest')
    subs2use = flist(contains(flist, 'rest'));
elseif strcmpi(subsessions2plot, 'all')
    subs2use = flist;
else 
    for i = 1:length(subsessions2plot)
       subs2use(i) = flist(contains(flist, ['sub' num2str(subsessions2plot(i)) '_']));
    end
end

for session = 1:length(subs2use)    
    obhvIn = load(subs2use{session});
    obhvIn = obhvIn.obhv;
    obhv(session) = obhvIn;
end
clear obhvIn;

subsessionInds = [];
allBhv.markers = [];
allBhv.trialTypeData= [];
allBhv.s = []; 
nTotalTrials = [];

for ss = 1:length(subs2use)

    sessionID = obhv(ss).Info.blockname;

    % camera variables different across subjects
    if strcmpi(subject, 'Barney') || strcmpi(subject, 'Dexter')
        faceCameraStream = 'FT1dat';
        bhvSR = 70;

    elseif strcmpi(subject, 'Thor')
        faceCameraStream = 'FT2dat';
        bhvSR = 120;
    end
    
    pulseLength = length(obhv(ss).(faceCameraStream).time);
    nFramesExported = size(obhv(ss).(faceCameraStream).Markers,1);
    if nFramesExported > pulseLength
        disp(['Warning -- sessionID ' num2str(obhv(ss).Info.blockname) ' dropped frames'])
        disp(['nFramesExported from DLC = ' num2str(nFramesExported)]);
        disp(['pulse Length = ' num2str(pulseLength)]);
    end

    %%%%%%% extract DLC marker positions & labels from obhv
    %%%%%
    % 
    DLCin = obhv(ss).(faceCameraStream);
    DLCtimeBase = DLCin.time(1:nFramesExported);
    [DLCmarkers, markerLabels] = processDLCMarkers('DLCin', DLCin, 'markers2use', markers2use,'excludeTongue', true,'smoothSize',smoothSize,'plotFlag',plotDLCFlag);
    subsessionInds(ss+1) = size(DLCmarkers,1); % length of each subsession

    %%%% extract trials 
    obhvThis = takeSpacedBhvs(obhv(ss),minRest,minRestFlag); % remove trials w/ <500 ms rest prior move onset
    allTrials = ismember(cell2mat(obhvThis.evScore(:,5)),bhvs2plot);

    % get rid of chew trials inside faceExp sets
    if strcmpi(subject, 'Thor')  && cell2mat(obhvThis.evScore(1,6)) == 3
        bhv2lose = ~(cell2mat(obhvThis.evScore(:,5)) == 4); 
        allTrials = logical(allTrials.*bhv2lose);
    end

    % remove events whose trial duration is outside of framesExported
    fmin = tmin*bhvSR; % convert to frames
    fmax = tmax*bhvSR;
    
    allOnsets = cell2mat(obhvThis.evScore(:,1));
    theseTrials = (allOnsets >= (fmin )) & (allOnsets + fmax <= nFramesExported);
    theseTrials = logical(allTrials .* theseTrials);
    
    theseOnsets = allOnsets(theseTrials); % handscored index, in subsession
    
    if strcmpi(sessionID, 'Thor-171010-154538') 
        disp(['sessionID = ' sessionID '  -- adjusting onsets for video export']);
        theseOnsets = theseOnsets - 2000; 
        theseOnsets = theseOnsets((theseOnsets-fmin) >= 0);
    end
    
    trialTypes = cell2mat(obhvThis.evScore(theseTrials,5));
    nTotalTrials(ss,:) = hist(trialTypes,bhvs2plot);

    % extract trialData
    for tt = 1:length(theseOnsets)
        onsetMin = theseOnsets(tt) - fmin;
        onsetMax = theseOnsets(tt) + fmax;
        thisData = DLCmarkers(onsetMin:onsetMax,:);

        %%%% downsample the DLC tracking data
        DLCtimeBase = -tmin:(1/bhvSR):tmax;
        [thisData, downsampledTimeAxis] = rebinBhv(thisData, DLCtimeBase, newBinSize);

        thisTrialTypeData = repmat(trialTypes(tt), [size(thisData,1)],1);

        allBhv.markers = [allBhv.markers; thisData];
        allBhv.trialTypeData = [allBhv.trialTypeData; thisTrialTypeData];
        allBhv.s = [allBhv.s ; obhv(ss).Info.blockname];
        
    end
    
   

end

%%% perform decomposition 
%% Input is nTimepoints x nMarkers,
if centerFlag
     [coeff,score,latent, ~, explained] = pca(allBhv.markers); % centering (subtract column/cell means) done automatically by pca()
elseif centerNormalizeFlag
     [coeff,score,latent, ~, explained] = pca(allBhv.markers,'VariableWeights', 'variance'); 
else
    [coeff,score,latent, ~, explained] = pca(allBhv.markers,'Centered','off');
    disp('performing PCA on uncentered, non-normalized data!')
end


%% per expression type, per trial, extract X&Y positions from -tmin/tmax around onset 
%  bhvTraj(bhv).data = nTimepoints x nMarkers x nTrials 
for bhv = 1:length(bhvs2plot)

    theseTrials = (allBhv.trialTypeData) == bhvs2plot(bhv);
    trialDur = length(downsampledTimeAxis);
    nTrials = sum(nTotalTrials);

    bhvTraj(bhv).data = reshape(score(theseTrials,:),[trialDur, nTrials(bhv),size(DLCmarkers,2)]);
    
    bhvTraj(bhv).avg = squeeze(mean(bhvTraj(bhv).data,2));
    bhvTraj(bhv).sem = std(bhvTraj(bhv).data,[],2) ./ sqrt(nTrials);

end

%%% do the TSNE calculations

if runTSNE
    nPCs2use = find(cumsum(latent) ./sum(latent) >= 0.95,1,'first');
    tsneResults = tsne(score(:,1:nPCs2use));
    figure;
    h = gscatter(tsneResults(:,1), tsneResults(:,2),allBhv.trialTypeData,'rbg',[],12);
    legend('Threat','Lipsmack','Chew'); legend boxoff
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    title('Characterization of Three Facial Expressions');
end

%% per expression type, per trial, extract new projection of markers 
%  bhvTraj(bhv).data = nTimepoints x nMarkers x nTrials 
for bhv = 1:length(bhvs2plot)

    theseTrials = (allBhv.trialTypeData) == bhvs2plot(bhv);
    trialDur = length(downsampledTimeAxis);
    nTrials = sum(nTotalTrials);

    bhvTraj(bhv).data = reshape(score(theseTrials,:),[trialDur, nTrials(bhv),size(DLCmarkers,2)]);
    
    bhvTraj(bhv).avg = squeeze(mean(bhvTraj(bhv).data,2));
    bhvTraj(bhv).sem = std(bhvTraj(bhv).data,[],2) ./ sqrt(nTrials);

end

% original marker data, separated by bhv 
for bhv = 1:length(bhvs2plot)

    theseTrials = (allBhv.trialTypeData) == bhvs2plot(bhv);
    trialDur = length(downsampledTimeAxis);
    nTrials = sum(nTotalTrials);

    markerData2plot(bhv).data = reshape(allBhv.markers(theseTrials,:),[trialDur, nTrials(bhv),size(DLCmarkers,2)]);
    
    markerData2plot(bhv).avg = squeeze(mean(markerData2plot(bhv).data,2));
    markerData2plot(bhv).sem = squeeze((std(markerData2plot(bhv).data,[],2)) ./ sqrt(nTrials(bhv)));

end

meanPos = mean(allBhv.markers,1);

%% plotting
scalar = 3000; 
PCweights = coeff(:,1:10);
percentVar = explained; 
[customColormap] = continuousColorMap(length(downsampledTimeAxis));

taxis = downsampledTimeAxis;

clear DLCin; clear obhv; clear DLCmarkers;

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
    title(['PC' num2str(PC) ' : ' num2str(percentVar(PC)) '% variance'])

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

%%
fig = figure('Position', get(0, 'Screensize')); 
    for PC = 1:10
        subplot(2,5,PC)
        plot(downsampledTimeAxis,bhvTraj(1).avg(:,PC),'r','linew',6); hold on; ...
            plot(downsampledTimeAxis,bhvTraj(2).avg(:,PC),'b','linew',6)
            plot(downsampledTimeAxis,bhvTraj(3).avg(:,PC),'g','linew',6); hold off
        title(['PC=' num2str(PC)],'FontSize',20)
    end

legend('Thr','LS', 'Chew')
sgtitle(['Top 10 FaceMotion PCs; tmin = ' num2str(tmin) ' tmax = ' num2str(tmax)],'FontSize',24)


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

%%
function [bhvOut,downsampledTimeBaseVector] = rebinBhv(bhvIn, timeBaseOld, newBinsize)
    
downsampledTimeBaseVector = timeBaseOld(1):newBinsize:timeBaseOld(end);
bhvOut = zeros(length(downsampledTimeBaseVector), size(bhvIn,2));

for j = 1:size(bhvIn,2)
    bhvOut(:,j) = interp1(timeBaseOld, bhvIn(:,j), downsampledTimeBaseVector);
end

end

