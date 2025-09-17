%%%%%
%%%%% Generates and analyzes a simulated dataset 
%%%%% for your preserved v. non-preserved neural correlations across
%%%%% bhvs analysis


clear all ; close all 

bhvs2plot = [1 2 4];
defBhv = 2;
nGestures = length(bhvs2plot);
regions = {'regionA','regionB'};
modeA = 'independent'; modeB = 'independent';
N1 = 50; N2 = 25; % N1, neurons in preserved region, 'N2' neurons in non-preserved region
perturbA = 0.15; perturbB = 0.15;
nTimepoints = 75;

nClust = 5;


sim = simulateRankOrderDatasetFlexible('modeA', modeA, 'modeB', modeB, 'N1', N1, 'N2', N2, 'T', nTimepoints, 'Seed',123);

% Outputs:
%   sim.regionBlah 
%       .X{1..3}   : T x N1 time series (behavior 1..3)
%       .R{1..3}   : N1 x N1 Pearson correlation matrices
%       .Ktrue{1..3}: N1 x N1 true covariances used for sampling


% pull out N1 x N1 Pearson correlation matrices per gesture 
for gg = 1:nGestures
    signalCorrOverTime.regionA(gg).corr = sim.regionA.R{1,gg};
    signalCorrOverTime.regionB(gg).corr = sim.regionB.R{1,gg};
end

for rr = 1:length(regions)

    % cluster & sort, according to ONE ordering (unified sorting)
    corrMatDef = signalCorrOverTime.(regions{rr})(defBhv).corr - eye(size(signalCorrOverTime.(regions{rr})(defBhv).corr));
    dissimSignalDef = 1 - corrMatDef';

    % compute distances, cluster distances, get new sorted indices of
    % defined bhv
    ZsigDef = linkage(dissimSignalDef,'complete');
    TsigDef = cluster(ZsigDef,'maxclust',nClust);
    [~, defSortIndxSig] = sort(TsigDef);

    % stash the sorted, self-defining bhv correlation matrix
    signalCorrOverTime.(regions{rr})(defBhv).unifiedSort.data = dissimSignalDef(defSortIndxSig,defSortIndxSig);
    signalCorrOverTime.(regions{rr})(defBhv).selfSort.data = dissimSignalDef(defSortIndxSig,defSortIndxSig);
   
    % Build aligned correlation matrices (same neuron order for all behaviors)
    R1 = signalCorrOverTime.(regions{rr})(1).corr(defSortIndxSig, defSortIndxSig);
    R2 = signalCorrOverTime.(regions{rr})(2).corr(defSortIndxSig, defSortIndxSig);
    R3 = signalCorrOverTime.(regions{rr})(3).corr(defSortIndxSig, defSortIndxSig);

    % (Optional but recommended) zero the diagonal to avoid any stray NaNs
    R1(1:size(R1,1)+1:end) = 0;
    R2(1:size(R2,1)+1:end) = 0;
    R3(1:size(R3,1)+1:end) = 0;

    % plot for sanity
    figure; subplot 131; imagesc(R1'); colormap(bluewhitered); 
    subplot 132; imagesc(R2'); colormap(bluewhitered)
    subplot 133; imagesc(R3'); colormap(bluewhitered)
    sgtitle(['Region : ' regions{rr} '/ preserved ordering'])

    % Run your rank-order similarity function(s)
    % (Use whichever version you implemented: Spearman-only or Spearman+Kendall)
    out.(regions{rr}) = rankOrderSpearmanKendall3(R1,R2, R3, 'NPerm',1000,'NBoot',200,'Seed',777);

end

out.regionA.summary
out.regionB.summary

