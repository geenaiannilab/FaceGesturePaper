% written GI updated 211111
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
% 
% This script will then plot a bunch of things to visualize different
% aspects of this reduced neural state space. 


clear all; close all;
set(0,'defaultAxesFontSize',36)
set(0,'defaultAxesFontWeight','bold')

date = '210704';
subject ='Barney';
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis']);
%outdir = '/Users/geena/Dropbox/PhD/Thesis/figures/';

chls = [1:240];
region = {'All'};
subsessions2plot =[2:5 7:10];
bhvs2plot = [1 2 4];

% DO NOT CHANGE!
win = 0.03;  % in sec
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

if centerOnlyFlag
     [coeff,score,latent, ~] = pca(neuralInput); % centering (subtract column/cell means) done automatically by pca()
     input2plot = neuralInput - mean(neuralInput); %for plotting only, later
elseif centerNormalizeFlag
     [coeff,score,latent, ~] = pca(neuralInput,'VariableWeights', 'variance'); 
     input2plot = (neuralInput - mean(neuralInput)) ./ std(neuralInput,[],1);
else
    [coeff,score,latent, ~] = pca(neuralInput,'Centered','off');
    input2plot = neuralInput;
    disp('performing PCA on uncentered, non-normalized data!')
end

% explained variance 
explained = latent/sum(latent);

% plotting PCA results

% Trial-Average Response (PSTH); input to PCA
markerSize = 36;
figure; imagesc(neuralInput'); colormap(bluewhitered); 
ylabel('cell'); xlabel('time'); 
xticklabels = -1:0.1:1;
xticks = linspace(1, size(neuralInput,1), numel(xticklabels));
yticks = 1:size(neuralInput,2);
yticklabels = spikeLabels;
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', yticks, 'YTickLabel', yticklabels)
colorbar;
title(['Trial-Avg FR per Cell: bhvtype: ' num2str(bhvs2plot)])

figure; imagesc(input2plot'); colormap(bluewhitered); 
ylabel('cell'); xlabel('time'); 
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', yticks, 'YTickLabel', yticklabels)
colorbar;
title(['PREPROC Trial-Avg FR per Cell: bhvtype: ' num2str(bhvs2plot)])


% plot FR vs loadings to see if highFR neurons determine own PCs
figure;
for i = 1:length(bhvs2plot) 
    subplot(1,2,1); 
    plot(meanFRs(thrCells2use,i), coeff(:,1)', 'o'); hold on;
end
hold off
grid; xlabel('meanFR'); ylabel('loadings on PC1')

for i = 1:length(bhvs2plot)
    subplot(1,2,2); 
    plot(meanFRs(thrCells2use,i), coeff(:,2)', 'o'); hold on;
end
hold off
grid; xlabel('meanFR'); ylabel('loadings on PC2')
sgtitle('Mean FRs vs. loadings on PC1, PC2')

% loadings of each neuron on each PC
 fig = figure; imagesc(coeff(:,1:5));  colorbar; colormap(bluewhitered)
 ax = gca; 
 ax.FontSize = 10;
 ax.YTick = (1:size(neuralInput,2));
 ax.YTickLabels = spikeLabels;
 ax.XTick = (1:5)
 xlabel('PC','FontSize',18); ylabel('Cell ID','FontSize',18);
 title('PC Loadings are Distributed Across Cells','FontSize',22)
%saveas(fig, [outdir '/' subject '_' date '_' region  '_weights_neuralTrajs.png']);
%savefig(fig, [outdir '/' subject '_' date '_' region  '_weights_neuralTrajs.fig']);

% eigenvalues 
fig = figure; 
plot(cumsum(latent/sum(latent)), 'ok','MarkerFaceColor','k') 
xlabel('PC'); ylabel('% variance'); xlim([0 80]);
title('Cumulative Variance Explained');
%saveas(fig, ['~/Desktop/neuralTrajs_cumVar_tmin' num2str(tmin2take) '.png']);
%saveas(fig, [outdir '/' subject '_' date '_' region  '_elbow_neuralTrajs.png']);
%savefig(fig, [outdir '/' subject '_' date '_' region  '_elbow_neuralTrajs.fig']);

% plot per TimePt data projected onto first 3 PCs, labeled by timebin
colorAxis = repmat(taxis2take,numel(score(:,1)),1);
colorAxis = colorAxis(1,:); markerSize = 36;
fig = figure; 
start = 1;
for i = 1:length(bhvs2plot)
    stop = length(taxis2take)*i;
    scatter3(score(start:stop,1),score(start:stop,2),score(start:stop,3),markerSize, colorAxis,'Linewidth',4)
    c = colorbar;
    set(c, 'ylim', [edges(3) edges(end-2)])
    hold on
    start = start + length(taxis2take);
end
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
title(['Time evolution of Neural Trajectories during Facial Expressions']); 
hold off; 
grid on; 
%saveas(fig, ['~/Desktop/neuralTrajs_tmin' num2str(tmin2take) '.png']);
%saveas(fig, [outdir '/' subject '_' date '_' region  '_TimeNeuralTrajs.png']);
%savefig(fig, [outdir '/' subject '_' date '_' region  '_TimeNeuralTrajs.fig']);

%%
% perTimePt data projected onto PCs 1-3; labeled by bhvtype
fig  = figure; 
markerType = 'osd*.';
markerSize = 300;
[customColormap] = continuousColorMap(length(taxis2take));

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
title(['Facial Gesture Neural Trajectories']); 
%legend('Threat','Lipsmack','Chew')
grid on; axes = gca; axes.LineWidth = 3; axes.GridLineStyle = ':';
%saveas(fig, [outdir '/' subject '_' date '_' region  '_neuralTrajs.png']);
%savefig(fig, [outdir '/' subject '_' date '_' region  '_neuralTrajs.fig']);
%%

% visualize each PC, for each bhv 
%%
figure; 
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
