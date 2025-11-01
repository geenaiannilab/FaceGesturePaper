function [compiled, taxis_out, gestureNames_out, bhvs2plot_out, meta] = ...
    compileDLCdays(dayRoots, varargin)
% Compile multiple days of DLC marker trajectories into gesture-wise norms.
% Each day folder must contain a 'DLC' subfolder holding the marker file.
%
% Modes:
%   'norm'      (default): compute Euclidean norms across *all* markers present per day.
%   'intersect'          : if 'markerNames' exist, restrict to the intersection of names
%                          (both X and Y for each marker) before computing norms.
%
% Outputs:
%   compiled(b).norms : [time x trialsAll] Euclidean norm trajectories for gesture b
%   taxis_out         : common time axis (must match across days)
%   gestureNames_out  : cellstr gesture names
%   bhvs2plot_out     : indices of gestures to plot (if present)
%   meta.files        : file paths used per day
%   meta.trialsPerDay : [nBehaviors x nDays] trial counts

    p = inputParser;
    addParameter(p,'Filename','markerTrajectories_10ms.mat',@(s)ischar(s)||isstring(s));
    addParameter(p,'Mode','norm',@(s)any(strcmpi(s,{'norm','intersect'})));
    parse(p,varargin{:});
    fnamePreferred = p.Results.Filename;
    mode = lower(p.Results.Mode);

    assert(iscell(dayRoots) && ~isempty(dayRoots), 'Provide dayRoots as a non-empty cell array of day directories.');

    % Holders
    compiled = []; taxis_out = []; gestureNames_out = []; bhvs2plot_out = [];
    meta = struct('files', {{}}, 'trialsPerDay', []);

    % For 'intersect' mode
    haveNamesAllDays = true;
    allNameSets = {};  % per day markerNames (without _x/_y suffix handling left to you)

    % First pass: load & check headers, gather common info
    dayData = cell(numel(dayRoots),1);
    for d = 1:numel(dayRoots)
        dlcDir = fullfile(dayRoots{d}, 'DLC');
        fPreferred = fullfile(dlcDir, fnamePreferred);
        fFallback  = fullfile(dlcDir, 'markerTrajectories.mat');
        if exist(fPreferred,'file'), fpath = fPreferred;
        elseif exist(fFallback,'file'), fpath = fFallback;
        else, error('File not found in %s', dlcDir);
        end
        S = load(fpath);  % expects at least: bhvTraj, taxis, gestureNames

        for rv = ["bhvTraj","taxis","gestureNames"]
            assert(isfield(S,rv), 'Missing variable "%s" in %s', rv, fpath);
        end

        % Initialize or consistency-check
        if isempty(taxis_out)
            taxis_out        = S.taxis(:);
            gestureNames_out = S.gestureNames;
            if isfield(S,'bhvs2plot'), bhvs2plot_out = S.bhvs2plot; end
        else
            assert(isequal(S.gestureNames, gestureNames_out), 'gestureNames differ across days.');
            assert(numel(S.taxis)==numel(taxis_out) && all(S.taxis(:)==taxis_out(:)), 'taxis differs across days.');
            if isfield(S,'bhvs2plot') && ~isempty(bhvs2plot_out)
                assert(isequal(S.bhvs2plot, bhvs2plot_out), 'bhvs2plot differs across days.');
            elseif isfield(S,'bhvs2plot') && isempty(bhvs2plot_out)
                bhvs2plot_out = S.bhvs2plot;
            end
        end

        % Keep data & metadata
        dayData{d} = S;
        meta.files{end+1} = fpath; %#ok<AGROW>

        % For intersect mode, record names
        if strcmp(mode,'intersect')
            if isfield(S,'markerNames') && ~isempty(S.markerNames)
                allNameSets{end+1} = S.markerNames(:); %#ok<AGROW>
            else
                haveNamesAllDays = false;
            end
        end
    end

    % If intersect mode requested but names missing, fall back to norm
    if strcmp(mode,'intersect') && ~haveNamesAllDays
        warning('markerNames not found in all days; falling back to ''norm'' mode.');
        mode = 'norm';
    end

    % Build intersection index per day if needed
    idxByDay = [];
    if strcmp(mode,'intersect')
        % Strict intersection of marker names across all days
        common = allNameSets{1};
        for k = 2:numel(allNameSets), common = intersect(common, allNameSets{k}, 'stable'); end
        assert(~isempty(common), 'No common markers across days.');
        % For each day, build indices picking [marker_x, marker_y] pairs if that’s the naming scheme.
        % Here we assume markerNames enumerates components in order; we’ll just use "ismember".
        idxByDay = cell(numel(dayData),1);
        for d = 1:numel(dayData)
            names = dayData{d}.markerNames(:);
            [tf,loc] = ismember(common, names);
            idxByDay{d} = loc(tf);  % indices of common names in this day
        end
    end

    % Initialize compiled per-behavior
    nB = numel(gestureNames_out);
    compiled = repmat(struct('norms', []), nB, 1);
    meta.trialsPerDay = zeros(nB, numel(dayData));

    % Helper: NaN-safe Euclidean norm across components
    nanNorm = @(X,dim) sqrt(nansum(X.^2, dim));

    % Second pass: compute norms per day, per behavior, then concat trials
    for d = 1:numel(dayData)
        S = dayData{d};
        for b = 1:nB
            X = S.bhvTraj(b).data;   % [time x trials x comps]
            if strcmp(mode,'intersect')
                if ~isempty(idxByDay{d})
                    X = X(:,:,idxByDay{d});  % keep only common components
                end
            end
            [T,Tr,C] = size(X);
            % Compute norms across components for each trial
            norms_d = zeros(T, Tr);
            for j = 1:Tr
                % NaN-safe norm across the third dimension
                % (handles differing component counts and possible NaNs)
                norms_d(:,j) = nanNorm(squeeze(X(:,j,:)), 2);
            end
            % Concatenate trials across days
            if isempty(compiled(b).norms)
                compiled(b).norms = norms_d;
            else
                compiled(b).norms = [compiled(b).norms, norms_d]; %#ok<AGROW>
            end
            meta.trialsPerDay(b,d) = Tr;
        end
    end

    if isempty(bhvs2plot_out), bhvs2plot_out = 1:nB; end
end
