% written GI updated 220104
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


close all;
clear all;
set(0,'defaultAxesFontSize',20)

subject2analyze ={'combined'};
nIterations = 50;
popType = 'unbalanced';
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/Pseudopopulations/' popType '/']);

% DO NOT CHANGE!
win = 0.01; % in sec
tmin = 1;
tmax = 1;
taxis = -tmin:win:tmax;  

centerOnlyFlag = 0; % ONLY mean center input to PCA (dont divide by std)
centerNormalizeFlag = 1; % normalize cells/variables by their variance in addition to mean-centering (PCA on correlation matrix)    
threshlowFRFlag = 1; % 0/1;  remove low FR neurons 
threshlowFR = 0.5; %  sp/s threshold 
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
edges = -tmin:win:tmax;  % use 'edges' to rebin (before conv)

% get rid of the edges where convolution results in artifacts & only use
% those bins later 
tmin2take = tmin - (extraBins*win); 
tmax2take = tmax - (extraBins*win);
taxis2take = -tmin2take:win:tmax2take;
taxis2takeIdx = edges >= taxis2take(1) & edges <= taxis2take(end);
        
% 
originalTaxis = -1:0.001:1.2;

%%%%%%%%%%%%%%%%%%
for ss = 1:length(subject2analyze)

    thisSubj = subject2analyze{ss};

    % repeat for each pseudopopulation
    for iter = 1:nIterations
        
        % load the pseudopopulation
        pseudoPop = load([workdir 'twoSubjPseudopopulation_N' num2str(iter) '.mat']);

        % define some variables
        arrayList = pseudoPop.regionList;
        bhvs2plot = unique(pseudoPop.pseudoPopulation.(subject2analyze{:}).bhvTrials);
        trialLabels = pseudoPop.pseudoPopulation.(thisSubj).bhvTrials;

        for aa = 1:length(arrayList)

            nCells = size(pseudoPop.pseudoPopulation.(thisSubj).(arrayList{aa}),3);

            % dimensions are ntrials x nbins x ncells
            allCellsPSTH.(arrayList{aa}) = zeros(length(trialLabels),length(edges)-1,nCells);

            for trial = 1:size(allCellsPSTH.(arrayList{aa}),1) %per trial
                for cell = 1:size(allCellsPSTH.(arrayList{aa}),3) % per cell
                        
                    raster = pseudoPop.pseudoPopulation.(thisSubj).(arrayList{aa})(trial,:,cell);
                    rasterTimes = originalTaxis(find(raster));
                    newBinnedSpikes = histcounts(rasterTimes, edges); % bin w/ extras for discarding conv artifact

                    allCellsPSTH.(arrayList{aa})(trial, :, cell) = conv(newBinnedSpikes, gaussKern, 'same') ./win; %divide by win for FR
                end
            end
            
             % remove low FR neurons
            if threshlowFRFlag
                meanFRs = squeeze(mean(allCellsPSTH.(arrayList{aa}),[1 2]));
                allCellsPSTH.(arrayList{aa}) = allCellsPSTH.(arrayList{aa})(:,:,meanFRs' >= threshlowFR);
            end

            %  average FRs; this is the trial-averaged response PER cell
            for i = 1:length(bhvs2plot)
                theseTrials = (pseudoPop.pseudoPopulation.(thisSubj).bhvTrials == bhvs2plot(i));
                allCellsTrialAvgResponse(i).bhv = squeeze(mean(allCellsPSTH.(arrayList{aa})(theseTrials, taxis2takeIdx, :)));
                allCellsTrialAvgResponse(i).nReps = sum(pseudoPop.pseudoPopulation.(thisSubj).bhvTrials == bhvs2plot(i));
            end

            % concatente all averaged responses for PCA input 
            neuralInput = vertcat(allCellsTrialAvgResponse(:).bhv);
            nFeatures = size(neuralInput,2);

            % Trial-Averaged PCA 
            % Input is nBins x nCells
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
            
            figure(1); 
            figure(2); 
            figure(3); 
%             % Trial-Average Response (PSTH); input to PCA
%             figure; imagesc(input2plot'); colormap(bluewhitered);
%             ylabel('cell'); xlabel('time');
%             set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
%             colorbar;
%             title(['PREPROC Trial-Avg FR per Cell: bhvtype: ' num2str(bhvs2plot')])

            % loadings of each neuron on each PC
%             fig = figure; imagesc(coeff(:,1:10));  colorbar; colormap(bluewhitered)
%             ax = gca;
%             ax.FontSize = 10;
%             ax.YTick = (1:size(neuralInput,2));
%             xlabel('PC','FontSize',18);
%             title('PC Loadings are Distributed Across Cells','FontSize',22)

%             % eigenvalues
%             fig = figure;
%             plot(cumsum(latent/sum(latent)), 'ok','MarkerFaceColor','k')
%             xlabel('PC'); ylabel('% variance'); xlim([0 nFeatures+1]);
%             title('Cumulative Variance Explained');

            % plot per TimePt data projected onto first 3 PCs, labeled by timebin
            colorAxis = repmat(taxis,numel(score(:,1)),1);
            colorAxis = colorAxis(1,:); markerSize = 36;
            
            figure(1); 
            start = 1;
            for i = 1:length(bhvs2plot)
                stop = length(taxis)*i;
                scatter3(score(start:stop,1),score(start:stop,2),score(start:stop,3),markerSize, colorAxis,'Linewidth',4)
                c = colorbar;
                % set(c, 'ylim', [edges(3) edges(end-2)])
                hold on
                start = start + length(taxis);
            end
            xlabel('1st Principal Component')
            ylabel('2nd Principal Component')
            zlabel('3rd Principal Component')
            title(['Time evolution of Neural Trajectories during Facial Expressions: ' arrayList{aa} ' ' popType]);
            grid on;

            % perTimePt data projected onto PCs 1-3; labeled by bhvtype
            figure(2);
            colors = [1 0 0 0.2; 0 0 1 0.2; 0 1 0 0.2];
            start = 1;

            for i = 1:length(bhvs2plot)
                stop = length(taxis)*i;
                scatter3(score(start:stop,1),score(start:stop,2),score(start:stop,3), [], colors(i,1:3),'MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4,'LineWidth',2);
                hold on
                start = start + length(taxis);
            end
            xlabel('1st Principal Component')
            ylabel('2nd Principal Component')
            zlabel('3rd Principal Component')
            title(['Facial Expression Neural Trajectories ' arrayList{aa} ' ' popType]);
            legend('Threat','Lipsmack','Chew')
            grid on;

            % visualize each PC, for each bhv
            figure(3);
            for PC = 1:5
                start = 1;
                for i = 1:length(bhvs2plot)
                    stop = length(taxis)*i;
                    PrincipleAxis = score(start:stop,PC);
                    subplot(1,5,PC); plot(taxis,PrincipleAxis, 'Color', colors(i,:), 'LineWidth', 2); xlabel('s'); ylabel('a.u.'); title(['PC' num2str(PC)'])
                    hold on
                    start = start + length(taxis);
                end
            end
            sgtitle(['PCs ' arrayList{aa} ' ' popType])
        clearvars -except pseudoPop subject2analyze thisSubj nIterations workdir...
            win tmin tmax taxis taxis2takeIdx edges originalTaxis minRestFlag minRest ...
            centerOnlyFlag centerNormalizeFlag threshlowFRFlag threshlowFR...
            gaussKern ss iter aa arrayList bhvs2plot trialLabels popType;


        end % end arrayList 

    end % end iterations 
end % end subjects 



% plotting PCA results



% calculate & plot how far away in neural space starting positions of each trajectory are 
% start = 1;
% for bhv = 1:length(bhvs2plot)
%     
%     stop = length(taxis2take)*bhv;
%     neuralTraj(bhv).data = score(start:stop,:);
%     
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
% xlim([0 60]); ylim([0 2.5]); 
% xlabel('dim'); ylabel('au'); legend('ThrvChew','LSvChew','ThrvLS');
% title('Distance in starting points of Neural trajectories, AFO dimensionality')
% subtitle(['Tmin = ' num2str(tmin2take) ' Tmax = ' num2str(tmax2take)])
% %saveas(fig, ['~/Desktop/neuralTrajsDistance_tmin' num2str(tmin2take) '.png']);



    
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
