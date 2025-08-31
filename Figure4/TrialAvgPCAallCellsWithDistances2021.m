% written GI updated 250817
%
% PUPRPOSE: Calculate & plot reduced neural "state space"
%   where each point in reduced neural state-space represents a BIN in TIME
%   
% Input to PCA is allCellsTrialAvgResponse, which is nBins x nCells, where
% each entry is the smoothed FR at that timepoint for that cell 
%
% This script will calculate the 1) loadings of each neuron onto each PC 
%                                2) the projections of the data onto the
%                                PCs
%                                3) variance explained by each PC 
%                                4) distance between neural trajectories
%                                over time
% 
% This script will then plot a bunch of things to visualize different
% aspects of this reduced neural state space. 


clear all; close all;
set(0,'defaultAxesFontSize',36); set(0, 'DefaultLineLineWidth', 4);
set(0,'defaultAxesFontWeight','bold')

date = '210704';
subject ='Barney';
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis']);
saveFlag = true; 

regions = getChannel2CorticalRegionMapping(subject, 1);

chls = 1:240;
subsessions2plot = [2:5 7:10];
bhvs2plot = [1 2 4];

% DO NOT CHANGE!
win = 0.02;  % in sec
tmin = 1;
tmax = 1;
minRestFlag = 1;
minRest = abs(tmin); % in sec; minimal rest prior to move onset (trials to include)

centerOnlyFlag = 0; % ONLY mean center input to PCA (dont divide by std)
centerNormalizeFlag = 1; % normalize cells/variables by their variance in addition to mean-centering (PCA on correlation matrix)    
threshlowFRFlag = 1; % 0/1;  remove low FR neurons 
threshlowFR = 0.1; %  sp/s threshold 
%%% NOTE RE: PREPROCESSING INPUT
%%% pca() in matlab does mean-centering *by default*, but does not divide by stdev by default 

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

restSubs = subs2show(contains(subs2show,'rest'));

% concatenate subsession spike times 
datAllSubs = [];
    
for session = 1:length(subs2show)   
    
    dat = [];
    obhvIn = load(subs2show{session});
    obhvIn = obhvIn.obhv;

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

if strcmpi(subject, 'Thor')
    [allSpikesAllSubs, allTrials] = getRidOfChewTrials(allSpikesAllSubs, allTrials);
end

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

% concatente all averaged responses for PCA input 
neuralInput = vertcat(allCellsTrialAvgResponse(:).bhv);

% remove low FR neurons 
if threshlowFRFlag
    thrCells2use = ~sum(meanFRs >= threshlowFR,2) == 0;
    neuralInput = neuralInput(:,thrCells2use);
    spikeLabels = spikeLabels(thrCells2use,:);
    spikeLabels2plot = spikeLabels2plot(thrCells2use,:);
end
   
% double check you have enough cells
 if size(neuralInput,2) <= 5
    disp(['fewer than 5 cells, (' num2str(sum(cells2use)) ') quitting'])
    return
 end

% Trial-Averaged PCA 
% Input is nBins x nCells,
%   where each entry is the trial-averaged smoothed FR within a time-bin
% Each point in state-space is a TIMEPOINT

for rr = 1:length(regions)

    % chls in this array
    chls2take = ismember(spikeLabels2plot(:,1), regions{rr}.channels);

    pcaInput = neuralInput(:,chls2take);

    if centerOnlyFlag
        [coeff,score,latent, ~] = pca(pcaInput); % centering (subtract column/cell means) done automatically by pca()
        input2plot = pcaInput - mean(pcaInput); %for plotting only, later
    elseif centerNormalizeFlag
        [coeff,score,latent, ~] = pca(pcaInput,'VariableWeights', 'variance');
        input2plot = (pcaInput - mean(pcaInput)) ./ std(pcaInput,[],1);
    else
        [coeff,score,latent, ~] = pca(pcaInput,'Centered','off');
        input2plot = pcaInput;
        disp('performing PCA on uncentered, non-normalized data!')
    end

    % store in structure

    allDataOut.(regions{rr}.label).coeff = coeff;
    allDataOut.(regions{rr}.label).score = score;
    allDataOut.(regions{rr}.label).latent = latent;
    allDataOut.(regions{rr}.label).explained = latent/sum(latent);
    allDataOut.(regions{rr}.label).nDimGreater90 = find(round(cumsum(allDataOut.(regions{rr}.label).explained),2) >= 0.9, 1);

    % pull out per bhv neural Traj
    start = 1;
    for bhv = 1:length(bhvs2plot)
        stop = length(taxis2take)*bhv;
        allDataOut.(regions{rr}.label).bhv(bhv).neuralTraj = score(start:stop,:);
        start = start + length(taxis2take);
    end

end



% calculate & plot how far away in neural space starting positions of each trajectory are

for rr = 1:length(regions)

    maxDim = allDataOut.(regions{rr}.label).nDimGreater90;

    % position coordinates
    ThrCoords = allDataOut.(regions{rr}.label).bhv(1).neuralTraj(:,1:maxDim);
    LSCoords = allDataOut.(regions{rr}.label).bhv(2).neuralTraj(:,1:maxDim);
    ChewCoords = allDataOut.(regions{rr}.label).bhv(3).neuralTraj(:,1:maxDim);

    % euclidean distance as a function of time, from dim =1 to max dimensions
    for dd = 1:maxDim

        % nNeuralDimensional distance vector, per bhv
        sosThrChDistance = 0; sosThrLSDistance = 0; sosLSChDistance = 0 ;

        for dim = 1:dd
            % nNeuralDimensional distance vector, per bhv
            sosThrChDistance = sosThrChDistance + ((ThrCoords(:,dim) - ChewCoords(:,dim)) .^2);
            sosThrLSDistance = sosThrLSDistance + ((ThrCoords(:,dim) - LSCoords(:,dim)) .^2);
            sosLSChDistance = sosLSChDistance + ((LSCoords(:,dim) - ChewCoords(:,dim)) .^2);
        end

        % nNeuralDimensional distance vector, per bhv
        allDataOut.(regions{rr}.label).distanceThrVCh(:,dd) = sqrt(sosThrChDistance);
        allDataOut.(regions{rr}.label).distanceThrVLS(:,dd) = sqrt(sosThrLSDistance);
        allDataOut.(regions{rr}.label).distanceLSVCh(:,dd) = sqrt(sosLSChDistance);
        allDataOut.(regions{rr}.label).distanceAverage(:,dd) = mean([allDataOut.(regions{rr}.label).distanceThrVCh(:,dd), allDataOut.(regions{rr}.label).distanceThrVLS(:,dd),allDataOut.(regions{rr}.label).distanceLSVCh(:,dd)],2);

    end % end max dimensionality

    if saveFlag 
        save([workdir '/neuralTraj_' regions{rr}.label '.mat'], 'allDataOut',...
            'date','subject','subsessions2plot','bhvs2plot','win','tmin','tmax','taxis2take');
    end
end


%%%%%%%%%%
%%%%%%%% plotting PCA results
%%
for rr = 1:length(regions)
    
    latent = allDataOut.(regions{rr}.label).latent;
    score = allDataOut.(regions{rr}.label).score;
    maxDim = allDataOut.(regions{rr}.label).nDimGreater90;

    % eigenvalues
    fig = figure;
    plot(cumsum(latent/sum(latent)), 'ok','MarkerFaceColor','k')
    xlabel('PC'); ylabel('% variance'); xlim([0 80]);
    title(['Region ' regions{rr}.label 'Cumulative Variance Explained']);

    % plot per TimePt data projected onto first 3 PCs, labeled by timebin
%     colorAxis = repmat(taxis2take,numel(score(:,1)),1);
%     colorAxis = colorAxis(1,:); markerSize = 36;
%     fig = figure;
%     start = 1;
%     for i = 1:length(bhvs2plot)
%         stop = length(taxis2take)*i;
%         scatter3(score(start:stop,1),score(start:stop,2),score(start:stop,3),markerSize, colorAxis,'Linewidth',4)
%         c = colorbar;
%         set(c, 'ylim', [edges(3) edges(end-2)])
%         hold on
%         start = start + length(taxis2take);
%     end
%     xlabel('1st Principal Component')
%     ylabel('2nd Principal Component')
%     zlabel('3rd Principal Component')
%     title(['Region :' regions{rr}.label  'Time evolution of Neural Trajectories during Facial Expressions']);
%     hold off;
%     grid on;
%     %saveas(fig, ['~/Desktop/neuralTrajs_tmin' num2str(tmin2take) '.png']);
%     %saveas(fig, [outdir '/' subject '_' date '_' region  '_TimeNeuralTrajs.png']);
%     %savefig(fig, [outdir '/' subject '_' date '_' region  '_TimeNeuralTrajs.fig']);

    %%
    % perTimePt data projected onto PCs 1-3; labeled by bhvtype
   
    markerType = 'osd*.';
    markerSize = 300;
    [customColormap] = continuousColorMap(length(taxis2take));
    
    f1 = figure; set(f1, 'Position',[476 92 1022 774]);%axes(haPlot(2)); 
    %haPlot = tight_subplot(1,3,[.01 .03],[.1 .01],[.01 .01]);
    %haPlot = haPlot';    
    %axes(haPlot(1)); 
    
    start = 1;

    for i = 1:length(bhvs2plot)
        stop = length(taxis2take)*i;
        scatter3(score(start:stop,1),score(start:stop,2),score(start:stop,3),markerSize, colormap(customColormap(:,:,i)),'filled',markerType(i))
        hold on
        start = start + length(taxis2take);
    end
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    title(['Region: ' regions{rr}.label ' Facial Gesture Neural Trajectories']);
    %legend('Threat','Lipsmack','Chew')
    grid on; axes = gca; axes.LineWidth = 3; axes.GridLineStyle = ':';
    exportgraphics(gcf, ['~/Desktop/NeuralTrajs_' regions{rr}.label '.pdf'], 'ContentType', 'vector');
    savefig(gcf, ['~/Desktop/NeuralTrajs_' regions{rr}.label '.fig']);
    
    f2 = figure; set(f2, 'Position',[476 92 1022 774]);%axes(haPlot(2)); 

    plot(taxis2take,allDataOut.(regions{rr}.label).distanceThrVCh(:,3),'Color',[0.3010 0.7450 0.9330]); hold on;
    plot(taxis2take,allDataOut.(regions{rr}.label).distanceThrVLS(:,3),'Color',[0.4660 0.6740 0.1880])
    plot(taxis2take,allDataOut.(regions{rr}.label).distanceLSVCh(:,3),'Color',[0.6350 0.0780 0.1840]	)
    plot(taxis2take, allDataOut.(regions{rr}.label).distanceAverage(:,3),'k','linew',6);
    xlabel('Time'), ylabel('Distance, au')
    title(['Region: ' regions{rr}.label ', Distance between Trajs in 3D'])
   
    f3 = figure; set(f3, 'Position',[476 92 1022 774]);%axes(haPlot(2)); 
    plot(taxis2take,allDataOut.(regions{rr}.label).distanceThrVCh(:,maxDim),'Color',[0.3010 0.7450 0.9330]); hold on;
    plot(taxis2take,allDataOut.(regions{rr}.label).distanceThrVLS(:,maxDim),'Color',[0.4660 0.6740 0.1880])
    plot(taxis2take,allDataOut.(regions{rr}.label).distanceLSVCh(:,maxDim),'Color',[0.6350 0.0780 0.1840]	)
    plot(taxis2take,allDataOut.(regions{rr}.label).distanceAverage(:,maxDim),'k','linew',6);
    xlabel('Time'), ylabel('Distance, au')
    title(['Region: ' regions{rr}.label ', Distance between Trajs, maxDim'])
    legend('Thr-Ch','Thr-LS','LS-Ch','Avg')
    exportgraphics(gcf, ['~/Desktop/DistanceBtwnNeuralTrajs_' regions{rr}.label '.pdf'], 'ContentType', 'vector');

    % visualize each PC, for each bhv
    %%
    figure; title(['Region : ' regions{rr}.label])
    colors = 'rbg';
    for PC = 1:3
        start = 1;
        for i = 1:length(bhvs2plot)
            stop = length(taxis2take)*i;
            PrincipleAxis = score(start:stop,PC);
            subplot(1,3,PC);
            plot(taxis2take,PrincipleAxis, colors(i), 'LineWidth', 6); xlabel('Time (sec)'); ylabel('a.u.'); title(['PC' num2str(PC)'])
            hold on
            start = start + length(taxis2take);
        end
    end


    %% TO make movie of time evolving neural trajectories...
    %
    % bhvTraj = reshape(score,[length(taxis2take) 3 size(coeff,1)]);
    %
    % fig = figure('units','normalized','outerposition',[0 0 1 1]);
    % set(gcf,'color','w');
    % colors = 'rbgk';
    % xlim([-5.8810   10.0000]);
    % ylim([-10.0000    5.4676]);
    % zlim([ -10.0000    5.3041]);
    % view(3); hold on
    % xlabel('1st Principal Component')
    % ylabel('2nd Principal Component')
    % zlabel('3rd Principal Component')
    % title(['Facial Expression Neural Trajectories']);
    % %legend('Threat','Lipsmack','Chew')
    % grid on;
    % %
    % % loops = size(bhvTraj,1);
    % % F(loops) = struct('cdata',[],'colormap',[]);
    % %
    % % % bhvTraj is nTimebins x nBhvs x nPCs
    % % for ii = 1:loops
    % %     for bhv = 1:length(bhvs2plot)
    % %         scatter3(bhvTraj(ii,bhv,1),bhvTraj(ii,bhv,2),bhvTraj(ii,bhv,3),markerSize, colors(bhv), 'LineWidth',4)
    % %         hold on;
    % %     end
    % %
    % %     pause(0.001);
    % %     F(ii) = getframe(gcf);
    % %
    % % end
    % %
    % % myVideo = VideoWriter(['~/Desktop/' subject '-' date '-NeuralTrajs.avi']);
    % % myVideo.FrameRate = 30;  % reduced from default of 30
    % % open(myVideo);
    % % writeVideo(myVideo, F);
    % % close(myVideo)
    %
    % %
    % % % calculate & plot how far away in neural space starting positions of each trajectory are
    % % start = 1;
    % % for bhv = 1:length(bhvs2plot)
    % %
    % %     stop = length(taxis2take)*bhv;
    % %     neuralTraj(bhv).data = score(start:stop,:);
    % %
    %     start = start + length(taxis2take);
    % end
    %
    %
    % for dim = 1:size(score,2)
    %     ThrCoords = neuralTraj(1).data(1,1:dim);
    %     LSCoords = neuralTraj(2).data(1,1:dim);
    %     ChewCoords = neuralTraj(3).data(1,1:dim);
    %
    %     dThrvChew(dim) = sqrt(sum((ThrCoords-ChewCoords).^2));
    %     dLSvChew(dim) = sqrt(sum((LSCoords-ChewCoords).^2));
    %     dThrvLS(dim) = sqrt(sum((ThrCoords-LSCoords).^2));
    % end
    % fig = figure; plot(1:dim, dThrvChew,'linew',4); hold on; plot(1:dim, dLSvChew,'linew',4); plot(1:dim,dThrvLS,'linew',4); hold off
    % %xlim([0 60]); ylim([0 2.5]);
    % xlabel('dim'); ylabel('au'); legend('ThrvChew','LSvChew','ThrvLS');
    % title('Distance in starting points of Neural trajectories, AFO dimensionality')
    % subtitle(['Tmin = ' num2str(tmin2take) ' Tmax = ' num2str(tmax2take)])
    % %saveas(fig, ['~/Desktop/neuralTrajsDistance_tmin' num2str(tmin2take) '.png']);
    %
    %
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
