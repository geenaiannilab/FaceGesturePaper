close all;
clear all;
set(0,'defaultAxesFontSize',36);
set(0,'defaultAxesFontWeight','bold');

date = '171027';
subject ='Thor';
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis']);
outDir = '~/Dropbox/subFigurePDFs/';
chls = 129:192; % DO NOT subselect at this stage, ANOVA was run on all 
subsessions2plot = [1 3 5 6 7];
bhvs2plot = [1 2 4];

colors = 'rbmgc';
saveFlag = 1;
frThresh = 0.1; %in hz

win = 0.001; % in sec
tmin = 1;
tmax = 1;

minRestFlag = 1;
minRest = abs(tmin); % in sec; minimal rest prior to move onset (trials to include)

normalizeFlag = 1; % by # trials per behavior

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate & then downsample gaussKern for convolution w/ spikes binned at "win" resolution
smParams.smoothSize = 50; %std in milisec
smParams.filterType = 1; % 1= gauss, 2= caus half-gauss, 3=box
smWin = genSmWin(smParams);
gaussKern = dwnsampleKernel(smWin,win);

extraBins = length(gaussKern); % extra bins to discard due to conv edges 
tmin = tmin +(extraBins*win); 
tmax = tmax + (extraBins*win);
edges = -tmin:win:tmax;  

% get rid of the edges where convolution results in artifacts & only use
% those bins later 
tmin2take = tmin - (extraBins*win); 
tmax2take = tmax - (extraBins*win);
taxis2take = -tmin2take:win:tmax2take;
bins2take = extraBins+1:length(taxis2take)+extraBins; 

% get relevant files 
flist = natsortfiles(dir([workdir '/bhvspikes_sub*.mat']));
flist = fullfile({flist.folder}, {flist.name}); 

% determine which type of subsession to plot
if strcmpi(subsessions2plot,'face expression')
    subs2show = flist(contains(flist, 'face expression'));
elseif strcmpi(subsessions2plot,'visual expression')
    subs2show = flist(contains(flist, 'visual expression'));
elseif strcmpi(subsessions2plot,'chew expression')
    subs2show = flist(contains(flist, 'chew expression'));
elseif strcmpi(subsessions2plot,'chew')
    subs2show = flist(contains(flist, 'chew'));
elseif strcmpi(subsessions2plot,'rest')
    subs2show = flist(contains(flist, 'rest'));
elseif strcmpi(subsessions2plot, 'all')
    subs2show = flist;
else 
    for i = 1:length(subsessions2plot)
       subs2show(i) = flist(contains(flist, ['sub' num2str(subsessions2plot(i)) '_']));
    end
end


% concatenate subsession spike times 
datAllSubs = [];
    
for session = 1:length(subs2show)   
    
    dat = [];
    obhvIn = load(subs2show{session});
    obhvIn = obhvIn.obhv;
    
    obhv = obhvIn; 
    obhv = takeSpacedBhvs(obhvIn,minRest,minRestFlag); % remove trials w/ <500 ms rest prior move onset
    clear obhvIn;
    
    counter = 1;        % output index 

    [evindx, ~] = find(cell2mat(obhv.evScore(:,5)) == bhvs2plot); % indices of desired behavior
    if ~isempty(evindx) 

            for ch = chls
                for unit = 1:length(obhv.spikes.el(ch).sp)
                    if ~isempty(obhv.spikes.el(ch).sp(unit).times)


                        onsets = cell2mat(obhv.evScore(evindx,1));
                        spTimes = obhv.spikes.el(ch).sp(unit).times;
                        
                        if strcmpi(subject,'Thor')
                            dat.binnedSpikes(:,:, counter) =   binSpikes( obhv.FT2dat.time(onsets), spTimes, win, tmin, tmax); %output is trials x bins
                        elseif strcmpi(subject,'Barney') || strcmpi(subject,'Dexter')
                            dat.binnedSpikes(:,:, counter) =   binSpikes( obhv.FT1dat.time(onsets), spTimes, win, tmin, tmax); %output is trials x bins
                        end
                        dat.name(counter,:) = string([num2str(ch) '-' num2str(unit)]);

                        counter = counter +1;

                    else % if no spiketimes for this unit, in this subsession, leave it empty

                        dat.binnedSpikes(counter,:,:) =  [];
                        dat.name(counter) = NaN;
                        counter = counter + 1;
                        
                    end
                end
            end
            
            % after loop thru all cells, add the immutable behavioral
            % data
            dat.frameON = cell2mat(obhv.evScore(evindx,1));
            dat.trialType = cell2mat(obhv.evScore(evindx,5));
            dat.condition = (obhv.evScore(evindx,6));
            dat.comments = obhv.evScore(evindx,8);

    else %if there's no relevant behavior in subsession, leave it empty
        dat.binnedSpikes = [];
        dat.name = [];
        dat.frameON = [];
        dat.trialType = [];
        dat.condition = [];
        dat.comments = [];
    end

    datAllSubs = [dat; datAllSubs]; %concatenate across subsessions
    
end
    
% concatenate all trials across all subsessions 
% dimensions are trials x bins x cells
allSpikesAllSubs.binnedSpikes = (vertcat(datAllSubs(:).binnedSpikes));   % divide by win to convert to FR
allSpikesAllSubs.name = datAllSubs.name;        

% behavioral / qualitative info (same for all cells) 
allTrials.frameON = vertcat(datAllSubs(:).frameON);
allTrials.trialType = vertcat(datAllSubs(:).trialType);
allTrials.condition = vertcat(datAllSubs(:).condition);
allTrials.comments = vertcat(datAllSubs(:).comments);

% dimensions are ntrials x nbins x ncells
allCellsPSTH = zeros(size(allSpikesAllSubs.binnedSpikes));

for trial = 1:size(allSpikesAllSubs.binnedSpikes,1) %per trial 
    for cell = 1:size(allSpikesAllSubs.binnedSpikes,3) % per cell 
        allCellsPSTH(trial, :, cell) = conv(allSpikesAllSubs.binnedSpikes(trial, :, cell), gaussKern, 'same') ./win; %divide by win for FR
    end
end


%  average FRs; this is the trial-averaged response PER cell 
%  allCellsTrialAvgResponse.bhv is nTimepoints x nCells (averaged over
%       nReps)
for i = 1:length(bhvs2plot)
    allCellsTrialAvgResponse(i).bhv = squeeze(mean(allCellsPSTH((allTrials.trialType == bhvs2plot(i)),bins2take,:)));
    allCellsTrialAvgResponse(i).nReps = sum(allTrials.trialType == bhvs2plot(i));
end

% mean FRs per cell, per bhv 
for bhv = 1:length(bhvs2plot)
    meanFRs(:,bhv) = mean(allCellsTrialAvgResponse(bhv).bhv,1);
end

% convert spike-labels to strings for plotting later
 spikeLabels = (allSpikesAllSubs(:).name);
for i= 1:length(spikeLabels) 
    spikeLabels2plot(i,:) = str2double(strsplit(spikeLabels(i),'-'));
end 

% load ANOVA results (-500 to 0 & 0 to 500+ ms)
%facePrefANOVA = load([workdir '\' date '_' subject '_facePrefAnova.mat']);

for unit =1:size(allCellsPSTH,3) % per cell 
    if sum(meanFRs(unit,:),2) <= frThresh
        continue
    end
    
    % get faceExp preference index 
    disp(['meanFRs= ' num2str(meanFRs(unit,:))]);
    [faceExpPref,depthPref] = calcFacePrefIndex(meanFRs(unit,:), length(bhvs2plot));
    fig1 = figure('units','normalized','outerposition',[0.3783    0.0835    0.6217    0.7974]);
    
     for bhv = 1: length(bhvs2plot) % per bhv
         
         thisBhvTrials = allTrials.trialType == bhvs2plot(bhv);
         
        %  fig1 = shadedErrorBar(taxis2take,(allCellsTrialAvgResponse(bhv).bhv(:,unit)),...
         %    {@(m) (mean(allCellsPSTH(thisBhvTrials,bins2take,unit))),...
          %   @(x) (std(allCellsPSTH(thisBhvTrials,bins2take,unit))) ./ sqrt(length(allCellsPSTH(thisBhvTrials,bins2take,unit)))},'lineprops',{'Color',colors(bhvs2plot(bhv))});
          fig1 = shadedErrorBar(taxis2take,(allCellsTrialAvgResponse(bhv).bhv(:,unit)),...
             {@(m) (mean(allCellsPSTH(thisBhvTrials,bins2take,unit))),...
@(x) (std(allCellsPSTH(thisBhvTrials,bins2take,unit))) ./ sqrt(size(allCellsPSTH(thisBhvTrials,bins2take,unit),1))},'lineprops',{'Color',colors(bhvs2plot(bhv)),'linew',8});
        
          hold on 
         set(gca,'box','off','tickdir','out');
         grid off
     end
     xline(0,'--','k','linew',5)
     hold off

     axes=gca;
     axes.FontWeight = 'bold'; Taxes.FontSize = 34;
     xlabel('Time, s');
     ylabel('Firing Rate, Hz');
     %legend({'Thr','LS','Chew'},'FontSize',17);
    
     title([subject(1) '-' date '-' spikeLabels{unit}]);
     subtitle(['FEPI = ' num2str(round(faceExpPref,1))],'FontSize',26);

    pause; 
    prompt = "Save? Y/N [Y]: ";
    txt = input(prompt,"s");
    if strcmpi(txt,'y')
        savefig(gcf, [outDir '/' subject '_' date '_PSTH_unit' spikeLabels{unit} '.fig'])
        %saveas(gcf, [outDir '/' subject '_' date '_PSTH_unit' spikeLabels{unit} '.png'])
        plotFig = gcf;
        exportgraphics(plotFig, [outDir '/' subject '_' date '_PSTH_unit' spikeLabels{unit} '.pdf'], 'ContentType', 'vector');
    else
        close all; 
    end
   
end


function binnedSpikes = binSpikes(evtime, sptimes, win, tmin, tmax)
%%%%%%%%%
ax = -tmin:win:tmax;
binnedSpikes = zeros(length(evtime),length(ax)-1); %initialize

for i=1:length(evtime)
    curT = evtime(i);
    spt = sptimes(find(sptimes >curT-tmin & sptimes < curT + tmax))-curT;%find spikes times after tmin & before tmax
                                                                         % subtracting curT converts spiketimes to -tmin:tmax scale 
                                                                         % (ie, first possible spktime is -tmin)
                                                                          
    binnedSpikes(i,:) = histcounts(spt, ax); %returns integer of # of spikes occurring in that bin


end
end
