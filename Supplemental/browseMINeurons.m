function browseMINeurons(binnedSpikes, labels, cellLabels, MI_corr, MI_raw, MI_shuff, pvals, binSize, taxis)
% Interactive browser: per-cell MI + smoothed time-resolved PSTH
%
% Inputs:
%   binnedSpikes : [nTrials x nTimeBins x nCells] binary or counts
%   labels       : [nTrials x 1] gesture identity (categorical or int)
%   MI_corr, MI_raw, MI_shuff, pvals : per-cell MI stats
%   binSize      : bin width in seconds (e.g. 0.05 for 50 ms bins)

    nCells = size(binnedSpikes,3);
    nGestures = numel(unique(labels));
    nTimeBins = size(binnedSpikes,2);
    timeAxis = taxis;

    % Gaussian kernel for smoothing (50 ms std)
    sigma = 0.05; % 50 ms
    winSize = round(5*sigma/binSize); % truncate at ±5σ
    gk = exp(-(-winSize:winSize).^2/(2*(sigma/binSize)^2));
    gk = gk / sum(gk);

    % Define colors explicitly
    gestureColors = [1 0 0; 0 0 1; 0 1 0]; % red, blue, green

    % Figure
    f = figure('Name','Neuron Browser','NumberTitle','off',...
               'KeyPressFcn',@keyPress);
    data.current = 1;
    data.nCells = nCells;
    data.binnedSpikes = binnedSpikes;
    data.labels = labels;
    data.MI_corr = MI_corr;
    data.MI_raw = MI_raw;
    data.MI_shuff = MI_shuff;
    data.pvals = pvals;
    data.nGestures = nGestures;
    data.timeAxis = timeAxis;
    data.gk = gk;
    data.binSize = binSize;
    data.gestureColors = gestureColors;
    data.cellLabels = cellLabels; 
    guidata(f,data);

    updatePlot(f);

    % --- Nested callbacks ---
    function keyPress(src,event)
        d = guidata(src);
        switch event.Key
            case {'rightarrow','downarrow','space'}
                d.current = min(d.current+1, d.nCells);
            case {'leftarrow','uparrow'}
                d.current = max(d.current-1, 1);
        end
        guidata(src,d);
        updatePlot(src);
    end

    function updatePlot(fig)
        d = guidata(fig);
        c = d.current;
        clf(fig);

        % Subplot 1: MI
        subplot(1,2,1); hold on;
        bar(1, d.MI_corr(c), 'FaceColor',[0.2 0.6 0.8]);
        errorbar(1, d.MI_raw(c), std(d.MI_shuff(:,c)), 'k.', 'LineWidth',1.5);
        plot(1, mean(d.MI_shuff(:,c)), 'rx', 'MarkerSize',10, 'LineWidth',2);
        if d.pvals(c) < 0.05
            text(1, d.MI_corr(c)+0.02, '*','HorizontalAlignment','center','FontSize',34);
        end
        xlim([0.5 1.5]);
        set(gca,'XTick',1,'XTickLabel',{sprintf('Neuron %d',c)});
        ylabel('MI (bits)');
        title('Mutual Information');
        legend('MI-c', 'MI-raw', 'mean, shuffle dist')

        % Subplot 2: time-resolved PSTH
        subplot(1,2,2); hold on;

        % Map labels to contiguous indices
        uniqueGestures = unique(d.labels);
        for g = 1:d.nGestures
            trials = d.labels == uniqueGestures(g);

            % average across trials for this gesture
            meanFR = squeeze(mean(d.binnedSpikes(trials,:,c),1)) / d.binSize; % Hz
            semFR  = squeeze(std(d.binnedSpikes(trials,:,c),[],1)) / sqrt(sum(trials)) / d.binSize;

            % smooth with Gaussian kernel
            meanFR_smooth = conv(meanFR, d.gk, 'same');
            semFR_smooth  = conv(semFR, d.gk, 'same');

            fill([d.timeAxis fliplr(d.timeAxis)], ...
                [meanFR_smooth+semFR_smooth fliplr(meanFR_smooth-semFR_smooth)], ...
                d.gestureColors(g,:), 'FaceAlpha',0.2, 'EdgeColor','none');
            plot(d.timeAxis, meanFR_smooth, 'Color',d.gestureColors(g,:), 'LineWidth',2);

            %shadedErrorBar(d.timeAxis,meanFR_smooth, 2*semFR_smooth,'lineprops',{'Color',d.gestureColors(g,:),'linew',8});

        
        end
        xlabel('Time (s)'); ylabel('Firing rate (Hz)');
        title(' PSTH ');
        legend(arrayfun(@(x) sprintf('Gesture %d',x), uniqueGestures,'UniformOutput',false));
        sgtitle(['Neuron  :', num2str(d.cellLabels(c,:))]);
    end
end
