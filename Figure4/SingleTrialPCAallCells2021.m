% written GI 
%%%     2021 NEEDS UPDATED HEADER COMMENTS !
%
% PUPRPOSE: Calculate & plot reduced neural "state space"
%   where each point in reduced neural state-space represents a TRIAL 
%   of a particular BEHAVIOR TYPE
%   
% Input to PCA is allCellsSpikeCounts, which is nCells x nTrials, 
%    where each entry is # of spikes that occured on that trial
%
% This script will calculate the 1) loadings of each neuron onto each PC 
%                                2) the projections of the data onto the
%                                PCs
%                                3) variance explained by each PC 
% 
% This script will then plot a bunch of things to visualize different
% aspects of this reduced neural state space. 


close all;
clear all;
set(0,'defaultAxesFontSize',24)

date = '210704';
subject ='Barney';
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis']);

chls = [33:128];
subsessions2plot = [2:5 7:10];%[8 9 10 11 12 13 14 ];%
bhvs2plot = [1 2 4];
minRestFlag = 1;
minRest = 0.5; % in sec; minimal rest prior to move onset (trials to include)

win = 0.02; % in sec
tmin = 1; tmax = 1; % 1.5 to -0.5; 1 to 0, 0.5 to 0.5
ax = -tmin:win:tmax;

sqrtFlag = 0;
centerOnlyFlag = 0; % ONLY mean center input to PCA (dont divide by std)
centerNormalizeFlag = 1; % normalize cells/variables by their variance in addition to mean-centering (PCA on correlation matrix)    
threshlowFRFlag = 1; % 0/1;  remove low FR neurons 
threshlowFR = 1; %  sp/s threshold 
%%% NOTE RE: PREPROCESSING INPUT
%%% pca() in matlab does mean-centering *by default*, but does not divide by stdev by default 

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

%%%%%%%%%%%%%
% sum spike counts over bins (now ntrials x ncells) 
allCellsPerTrialSpikeCount = squeeze(sum(allSpikesAllSubs.binnedSpikes,2));

% convert spike-labels to strings for plotting later
 spikeLabels = (allSpikesAllSubs(:).name);
for i= 1:length(spikeLabels) 
    spikeLabels2plot(i,:) = str2double(strsplit(spikeLabels(i),'-'));
end 


% calc Fano Factors 
bin.means = mean(allCellsPerTrialSpikeCount,1); % avg # of spikes/trial; NOT RATE
bin.SDs = std(allCellsPerTrialSpikeCount,[],1);
bin.meanFRs = mean(allCellsPerTrialSpikeCount,1)/(tmin+tmax); % now in sp/s

%plot FF
figure;
subplot(1,2,1); plot( bin.means, bin.SDs.^2,'o'); grid; line( [0 10], [0 10], 'linewidth',3,'color','k')
title('Fano factor of all Neurons')
tmp = (bin.meanFRs >= threshlowFR);
subplot(1,2,2); plot(bin.means(tmp), bin.SDs(tmp).^2, 'o');grid; line( [0 10], [0 10], 'linewidth',3,'color','k')
text(bin.means(tmp), bin.SDs(tmp).^2, spikeLabels(tmp))
title('Fano factor for Neurons w meanFR >1 hz')

% % remove low FR neurons 
if threshlowFRFlag
    allCellsPerTrialSpikeCount = allCellsPerTrialSpikeCount(:,(bin.meanFRs >= threshlowFR));
    spikeLabels = spikeLabels(bin.meanFRs >= threshlowFR,:);
end
   
% double check you have enough cells
 if size(allCellsPerTrialSpikeCount,2) <= 5
    disp(['fewer than 5 cells, (' num2str(sum(cells2use)) ') quitting'])
    return
 end


% figure; 
% colors = 'rbgck';
% for bhv = 1:length(bhvs2plot)
%     boxplot(allCellsPerTrialSpikeCount((allTrials.trialType==bhvs2plot(bhv)),:), 'Colors', colors(bhv), 'PlotStyle', 'compact');
%     ylabel('spikes per trial')
%     xticks(); xticklabels();
%     xticks(1:size(allCellsPerTrialSpikeCount,2));
%     xticklabels(spikeLabels);
%     xtickangle(45)
%     hold on
% end
% hold off;


% Single-Trial PCA 
% Input is nTrials x nCells , 
%    where each entry is # of spikes within a trial
%  Each point in state-space is a TRIAL 

if centerOnlyFlag
     [coeff,score,latent, ~] = pca(allCellsPerTrialSpikeCount); % centering (subtract column/cell means) done automatically by pca()
     neuralInput = allCellsPerTrialSpikeCount - mean(allCellsPerTrialSpikeCount); %for plotting only, later
elseif centerNormalizeFlag
     [coeff,score,latent, ~] = pca(allCellsPerTrialSpikeCount,'VariableWeights', 'variance'); 
     neuralInput = (allCellsPerTrialSpikeCount - mean(allCellsPerTrialSpikeCount)) ./ std(allCellsPerTrialSpikeCount,[],1);
elseif sqrtFlag 
    [coeff,score,latent, ~] = pca(sqrt(allCellsPerTrialSpikeCount)); 
    neuralInput = sqrt(allCellsPerTrialSpikeCount);
else
    [coeff,score,latent, ~] = pca(allCellsPerTrialSpikeCount,'Centered','off');
    neuralInput = allCellsPerTrialSpikeCount;
    disp('performing PCA on uncentered, non-normalized data!')
end

% explained variance 
explained = latent/sum(latent);

% plotting PCA results

% spikeCounts per cell per Trial (input to PCA)
markerSize = 36;
% figure; imagesc(allCellsPerTrialSpikeCount'); colormap(jet); ylabel('cell'); xlabel('trial#'); 
% colorbar;
% xticks(1:length(allTrials.trialType));
% xticklabels(allTrials.trialType);
% yticks(1:size(allCellsPerTrialSpikeCount,2));
% yticklabels(spikeLabels);
% title('spikes PER trial, per Cell ')
% 
% 
% figure; imagesc(neuralInput'); colormap(jet); ylabel('cell'); xlabel('trial#'); 
% colorbar;
% xticks(1:length(allTrials.trialType));
% xticklabels(allTrials.trialType);
% yticks(1:size(allCellsPerTrialSpikeCount,2));
% yticklabels(spikeLabels);
% title('spikes PER trial per cell, PREPROC INPUT ')
% 
% % % plot FR vs loadings to see if highFR neurons determine own PCs
% figure;
% subplot(1,2,1); plot(bin.meanFRs(bin.meanFRs >= threshlowFR), coeff(:,1), 'o'); grid; xlabel('meanFR'); ylabel('loadings on PC1')
% [rho1, pval1] = corr(bin.meanFRs(bin.meanFRs >= threshlowFR)', coeff(:,1));
% subplot(1,2,2); plot(bin.meanFRs(bin.meanFRs >= threshlowFR), coeff(:,2), 'o'); grid; xlabel('meanFR'); ylabel('loadings on PC2')
% [rho2, pval2] = corr(bin.meanFRs(bin.meanFRs >= threshlowFR)', coeff(:,2));
% sgtitle('Mean FRs vs. loadings on PC1, PC2')
% 

% loadings of each neuron on each PC
 figure; imagesc(coeff(:,1:5));  colorbar; colormap(bluewhitered)
 ax = gca; 
 ax.FontSize = 10;
 yticks(1:size(allCellsPerTrialSpikeCount,2));
 yticklabels(spikeLabels);
 xticks(1:5)
 xlabel('PC','FontSize',18); ylabel('Cell ID','FontSize',18);
 title('PC Loadings are Distributed Across Cells','FontSize',22)


% eigenvalues 
figure; 
plot(cumsum(latent/sum(latent)), 'ok','MarkerFaceColor','k') 
xlabel('PC'); ylabel('% variance'); xlim([0 80]);
title('Cumulative Variance Explained');

% perTrial data projected onto first 3 PCs 
%%
[~,~,id] = unique(allTrials.trialType);
colors = 'rbgck';
markers = 'osd*.';
plotFig = figure;
for idx = 1:length(unique(allTrials.trialType))'
    data = score(allTrials.trialType == bhvs2plot(idx),:);
    scat(idx) = scatter3(data(:,1), -data(:,2), data(:,3), 500, [colors(idx) markers(idx)], 'filled',...
        'MarkerFaceAlpha', 0.5,'MarkerEdgeColor',colors(idx),'linew',2);
    hold on;
end
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title(['Single Trials of Facial Gesture ']); 
legend('Threat','Lipsmack','Chew')
%%
% perTrial data projected onto 2PCs, labeled by bhvType 
% figure; 
% subplot(1,3,1); 
% gscatter(score(:,1),score(:,2),allTrials.trialType, colors(bhvs2plot), markers(bhvs2plot)); xlabel('PC1'); ylabel('PC2'); legend
% subplot(1,3,2);
% gscatter(score(:,2),score(:,3),allTrials.trialType, colors(bhvs2plot), markers(bhvs2plot)); xlabel('PC2'); ylabel('PC3'); legend
% subplot(1,3,3)
% gscatter(score(:,1),score(:,3),allTrials.trialType, colors(bhvs2plot), markers(bhvs2plot)); xlabel('PC1'); ylabel('PC3'); legend

% figure; 
% colorMap = [1,0,0;0,0,1;1,1,0;0,1,0;1,0,1];
% for PC = 1:3
%      PrincipleAxis = score(:,PC); 
%      subplot(1,3,PC); 
%      b(PC) = bar(PrincipleAxis);
%      b(PC).FaceColor = 'flat';
%      
%     for ii = 1:length(b(PC).CData)
%         b(PC).CData(ii,:) = colorMap(allTrials.trialType(ii),:); 
%     end
%     xlabel('trial#')
%     ylabel('weight on PC, a.u.')
%    
% end
% c = colorbar;
% c.Ticks = 0:length(bhvs2plot)-1;
% c.TickLabels = bhvs2plot;
% colormap(colorMap)
% sgtitle('Princple Components 1-3')
% sgtitle('Princple Components 1-3')
% 

% % biplots in 2D & 3D 
% 
% figure; h1 = biplot(coeff(:,1:3), 'scores', score(:,1:3), 'Color', 'b','varlabels', spikeLabels);
% title(['biplot PC1-3; bhvtype: ' num2str(bhvs2plot)]); 
% 
% figure; h2 =  biplot(coeff(:,1:2), 'scores', score(:,1:2), 'Color', 'b','varlabels', spikeLabels);
% title(['biplot PC1v2; bhvtype: ' num2str(bhvs2plot)]); 
% 
% figure; h3 = biplot(coeff(:,[1,3]), 'scores', score(:,[1,3]), 'Color', 'b','varlabels', spikeLabels);
% title(['biplot PC1v3; bhvtype: ' num2str(bhvs2plot)]); 
% 
% % labeling neurons by array & trials by bhvType
% %colorBiplotVectorsByArray(allCellsPerTrialSpikeCount', spikeLabels2plot, h1, h2, h3);
% %colorObservationsOnBiplot(allCellsPerTrialSpikeCount', allTrials.trialType, h1, h2, h3)
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

