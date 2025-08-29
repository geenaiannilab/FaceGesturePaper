
set(0,'defaultAxesFontSize', 24); % bc im blind
set(0,'defaultAxesFontWeight', 'bold'); 

dirs = {'Barney_210704', 'Barney_210706', 'Barney_210805','Thor_171005',...
    'Thor_171010', 'Thor_171027', 'Thor_171128'};

combined = summarize_rankorder_results_dirs(dirs, 'DayLabel','dir', ...
    'FilePattern','*rankSim*.mat', ...  % or '*.mat'
    'Recurse',false, ...
    'SaveCSV',true, ...
    'CSVPath','/Users/geena/Dropbox/PhD/SUAinfo/combined_rankorder_summary.csv');     

% plot spearman by region summary 
% plot_rankorder_by_region(combined, ...
%     'Metric','Spearman', ...
%     'RegionSrc','File', ...
%     'RegionPatterns', {'All','All'; 'S1','S1'; 'M1','M1'; 'PMv','PMv'; 'M3','M3'}, ...
%     'NBoot',5000, 'ShowPoints',true);

% Spearman summary
plot_rankorder_summary(combined, 'Metric','Spearman', 'Alpha',0.05, 'NBoot',5000, 'ShowPoints',true);

% Kendall summary (optional)
%plot_rankorder_summary(combined, 'Metric','Kendall', 'Alpha',0.05, 'NBoot',5000, 'ShowPoints',true);


function [combined] = summarize_rankorder_results_dirs(dirList, varargin)
% summarize_rankorder_results_dirs
% Walk a list of directories, load rank-order result tables, and print/aggregate.
%
% INPUT
%   dirList : string array or cellstr of directories (absolute or relative)
%
% Name-Value options (all optional):
%   'FilePattern' : '*.mat'             % which files to read in each dir
%   'Recurse'     : false               % search subfolders too (requires '**' support)
%   'SaveCSV'     : true                % write combined CSV
%   'CSVPath'     : 'combined_rankorder_summary.csv'
%   'DayLabel'    : 'dir'               % 'dir' (folder name) or 'file' (file base name)
%
% OUTPUT
%   combined : table with Day + all summary columns stacked across dirs/files
%
% Each result file must contain either:
%   - variable 'out' with field 'summary' (table), OR
%   - variable 'summary' directly (table)

% ------------- options -------------
p = inputParser;
addParameter(p,'FilePattern','*.mat',@(s)ischar(s)||isstring(s));
addParameter(p,'Recurse',false,@islogical);
addParameter(p,'SaveCSV',true,@islogical);
addParameter(p,'CSVPath','combined_rankorder_summary.csv',@(s)ischar(s)||isstring(s));
addParameter(p,'DayLabel','dir',@(s)any(strcmpi(string(s),["dir","file"])));
parse(p,varargin{:});
filePattern = string(p.Results.FilePattern);
recurse     = p.Results.Recurse;
saveCSV     = p.Results.SaveCSV;
csvPath     = string(p.Results.CSVPath);
daySource   = lower(string(p.Results.DayLabel));

% required columns
reqCols = {'Pair','Spearman_r','Spearman_p','Spearman_CI_low','Spearman_CI_high', ...
           'Kendall_tau','Kendall_p','Kendall_CI_low','Kendall_CI_high'};

% normalize dirList
if ischar(dirList) || isstring(dirList), dirList = cellstr(dirList); end
dirList = string(dirList(:));

allChunks = {};
fprintf('Scanning %d directories...\n\n', numel(dirList));

for d = 1:numel(dirList)
    dpath = ['/Users/geena/Dropbox/PhD/SUAInfo/' dirList(d) '/CorrMatricesV2/'];
    dpath = [dpath{:}];
    if ~isfolder(dpath)
        warning('Skipping: %s (not a folder)', dpath);
        continue;
    end

    % discover files
    if recurse
        S = dir(fullfile(dpath, '**', filePattern));
    else
        S = dir(fullfile(dpath, filePattern));
    end
    if isempty(S)
        warning('No files in %s matching %s', dpath, filePattern);
        continue;
    end

    % pretty day label (folder name)
    [~, dayLabelFromDir] = fileparts(char(dpath));

    fprintf('=== %s ===\n', dayLabelFromDir);

    for k = 1:numel(S)
        fn = fullfile(S(k).folder, S(k).name);
        data = load(fn);

        % pull summary table
        if isfield(data, 'out') && isfield(data.out, 'summary') && istable(data.out.summary)
            T = data.out.summary;
        elseif isfield(data, 'summary') && istable(data.summary)
            T = data.summary;
        else
            warning('  Skipping %s (no out.summary / summary table).', fn);
            continue;
        end

        % ensure required columns exist
        missing = setdiff(reqCols, T.Properties.VariableNames);
        if ~isempty(missing)
            warning('  Skipping %s (missing columns: %s)', fn, strjoin(missing,', '));
            continue;
        end

        % pull out dayLabel 
         dirLabel = string(dirList(d));  % e.g., 'Barney_210704'
        if daySource == "file"
            [~, baseName, ~] = fileparts(S(k).name);
            dayLabel = string(baseName);            % label by file basename
        else
            dayLabel = dirLabel;                    % label by dirList entry
        end

        % ---- print human-readable summary for this file/day ----
        fprintf('  File: %s\n', S(k).name);
        for i = 1:height(T)
            pr  = T.Spearman_r(i);
            psp = T.Spearman_p(i);
            cil = T.Spearman_CI_low(i);
            cih = T.Spearman_CI_high(i);

            kt  = T.Kendall_tau(i);
            pkp = T.Kendall_p(i);
            kil = T.Kendall_CI_low(i);
            kih = T.Kendall_CI_high(i);

            fprintf('    %-7s  ρ=% .3f (95%% CI [% .3f, % .3f]);  τ=% .3f (95%% CI [% .3f, % .3f])\n', ...
                T.Pair{i}, pr, cil, cih, kt, kil, kih);
            fprintf('             Mantel p(ρ)=%s   Mantel p(τ)=%s\n', pformat(psp), pformat(pkp));
        end
        fprintf('\n');

        % augment and stash
        Day   = repmat(dayLabel, height(T), 1);
        File  = repmat(string(S(k).name), height(T), 1);
        Chunk = addvars(T, Day, File, 'Before', 1);
        allChunks{end+1} = Chunk; %#ok<AGROW>
    end
end

if isempty(allChunks)
    warning('No valid result tables found. Returning empty table.');
    combined = table();
    return;
end

combined = vertcat(allChunks{:});
combined = movevars(combined, 'Day', 'Before', 1);

% save CSV
if saveCSV
    if ~exist(csvPath)
        fclose(fopen(csvPath, 'a'));  % "touch" the file: create if missing, then close
    end
    writetable(combined, csvPath);
    fprintf('Combined table saved to: %s\n', csvPath);
end
end

function plot_rankorder_summary(combined, varargin)
% plot_rankorder_summary  Collapse across days and plot summary with CIs & significance.
% Labels: R1='threat', R2='lipsmack', R3='chew'. Collapses over Day/File.

p = inputParser;
addParameter(p,'Metric','Spearman',@(s)ischar(s)||isstring(s));
addParameter(p,'Alpha',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'NBoot',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'ShowPoints',true,@islogical);
parse(p,varargin{:});
metric = lower(string(p.Results.Metric));
alpha  = p.Results.Alpha;
NBoot  = p.Results.NBoot;
showPts= p.Results.ShowPoints;

switch metric
    case "spearman"
        effCol = 'Spearman_r'; pCol = 'Spearman_p';
        ylab   = 'Rank-order similarity (Spearman \rho)';
    case "kendall"
        effCol = 'Kendall_tau'; pCol = 'Kendall_p';
        ylab   = 'Rank-order similarity (Kendall \tau)';
    otherwise
        error('Metric must be ''Spearman'' or ''Kendall''.');
end

needCols = {'Pair',effCol,pCol};
missing = setdiff(needCols, combined.Properties.VariableNames);
if ~isempty(missing), error('Combined table missing: %s', strjoin(missing,', ')); end

pairOrder   = {'R1-R2','R1-R3','R2-R3'};
pairDisplay = {'threat–lipsmack','threat–chew','lipsmack–chew'};
M = numel(pairOrder);

mu    = nan(M,1);
CI    = nan(2,M);
pComb = nan(M,1);
nPer  = nan(M,1);

for i = 1:M
    sel   = strcmp(combined.Pair, pairOrder{i});
    vals  = combined.(effCol)(sel);
    pvals = combined.(pCol)(sel);

    vals  = vals(isfinite(vals));
    pvals = pvals(isfinite(pvals));

    nPer(i) = numel(vals);
    if isempty(vals), continue; end

    mu(i) = median(vals, 'omitnan');

    if numel(vals)==1
        CI(:,i) = [vals; vals];
    else
        bs = nan(NBoot,1);
        for b = 1:NBoot
            ii = randi(numel(vals), numel(vals), 1);
            bs(b) = median(vals(ii), 'omitnan');
        end
        CI(:,i) = prctile(bs, [2.5 97.5])';
    end

    if ~isempty(pvals)
        pvals = max(pvals, 1e-12); % floor to avoid -Inf
        X = -2 * sum(log(pvals));
        df = 2 * numel(pvals);
        pComb(i) = 1 - chi2cdf(X, df);
    end
end

figure('Color','w'); hold on;
x = (1:M).';

% (optional) per-day points
if showPts
    gcol = 0.7*[1 1 1];
    for i = 1:M
        sel  = strcmp(combined.Pair, pairOrder{i});
        vals = combined.(effCol)(sel);
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end
        jitter = (rand(size(vals))-0.5)*0.15;
        plot(i + jitter, vals, '.', 'Color', gcol, 'MarkerSize', 24);
    end
end

% --- FIX: ensure column vectors for errorbar inputs ---
mu_col   = mu(:);                % Mx1
ci_low   = CI(1,:).';            % Mx1
ci_high  = CI(2,:).';            % Mx1
yneg     = mu_col - ci_low;      % Mx1
ypos     = ci_high - mu_col;     % Mx1

% Guard against negatives from NaNs/order issues (optional)
% yneg(yneg<0) = NaN; ypos(ypos<0) = NaN;

eh = errorbar(x, mu_col, yneg, ypos, 'o', 'LineWidth',3, 'CapSize',18, 'MarkerSize',24);

% significance stars
finiteCI = isfinite(ci_low) & isfinite(ci_high);
rngY = diff([min(ci_low(finiteCI)), max(ci_high(finiteCI))]);
if isempty(rngY) || ~isfinite(rngY) || rngY==0, rngY = 0.1; end
for i = 1:M
    if ~isfinite(pComb(i)), continue; end
    if pComb(i) < 0.001, stars='***';
    elseif pComb(i) < 0.01, stars='**';
    elseif pComb(i) < alpha, stars='*';
    else, stars='n.s.';
    end
    yTop = ci_high(i);
    if ~isfinite(yTop), yTop = mu_col(i); end
    text(x(i), yTop + 0.03*rngY, stars, 'HorizontalAlignment','center', 'FontSize', 16,'FontWeight','bold');
end

xlim([0.5 M+0.5]); ylim([-1 1])
xticks(x); xticklabels(pairDisplay);
ylabel(ylab);
title(sprintf('Across-day summary (%s) — Median with bootstrap 95%% CI; p via Fisher''s method', char(metric)));
grid on; box off;

annotation('textbox',[0.13 0.01 0.8 0.05], 'String', ...
    sprintf('Stars reflect combined Mantel p (Fisher). Points = per-day estimates (n per pair: %s).', ...
    strjoin(arrayfun(@(n)sprintf('%d',nPer(i)), (1:M), 'UniformOutput',false), ', ')), ...
    'EdgeColor','none','HorizontalAlignment','left','FontSize',14);
end

function plot_rankorder_by_region(combined, varargin)
% plot_rankorder_by_region  Make one summary plot per cortical region.
%
% Usage:
%   plot_rankorder_by_region(combined, ...
%       'Metric','Spearman', ...          % or 'Kendall'
%       'RegionSrc','File', ...           % which column to parse: 'File' or 'Day' or a custom col
%       'RegionPatterns', {'All','All'; 'S1','S1'; 'M1','M1'; 'PMv','PMv'; 'M3','M3'}, ...
%       'Alpha',0.05, 'NBoot',5000, 'ShowPoints',true)
%
% Notes:
%   - Collapses across days within each region.
%   - Pairs are displayed as: R1–R2 => threat–lipsmack, R1–R3 => threat–chew, R2–R3 => lipsmack–chew.
%   - Significance per pair is combined across days via Fisher’s method from per-day Mantel p-values.
%
p = inputParser;
addParameter(p,'Metric','Spearman',@(s)ischar(s)||isstring(s));
addParameter(p,'RegionSrc','File',@(s)ischar(s)||isstring(s));
addParameter(p,'RegionPatterns', {'All','All'; 'S1','S1'; 'M1','M1'; 'PMv','PMv'; 'M3','M3'});
addParameter(p,'Alpha',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'NBoot',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'ShowPoints',true,@islogical);
parse(p,varargin{:});
metric   = lower(string(p.Results.Metric));
regionSrc= char(p.Results.RegionSrc);
patMap   = p.Results.RegionPatterns;
alpha    = p.Results.Alpha;
NBoot    = p.Results.NBoot;
showPts  = p.Results.ShowPoints;

% ---- columns & labels ----
switch metric
    case "spearman", effCol='Spearman_r'; pCol='Spearman_p'; ylab='Rank-order similarity (Spearman \rho)';
    case "kendall",  effCol='Kendall_tau'; pCol='Kendall_p';  ylab='Rank-order similarity (Kendall \tau)';
    otherwise, error('Metric must be ''Spearman'' or ''Kendall''.');
end
needCols = {'Pair',effCol,pCol,regionSrc};
missing = setdiff(needCols, combined.Properties.VariableNames);
if ~isempty(missing), error('Combined table missing: %s', strjoin(missing,', ')); end

pairOrder   = {'R1-R2','R1-R3','R2-R3'};
pairDisplay = {'threat–lipsmack','threat–chew','lipsmack–chew'};
M = numel(pairOrder);

% ---- assign Region from pattern map (first match wins) ----
region = strings(height(combined),1);
srcVals = string(combined.(regionSrc));
for i = 1:size(patMap,1)
    pat   = string(patMap{i,1});          % pattern to search
    label = string(patMap{i,2});          % label to assign
    hit = contains(srcVals, pat, 'IgnoreCase', true);
    region(hit) = label;
end
% keep only rows that matched a known region
keep = region ~= "";
combined_region = combined(keep,:);
combined_region.Region = region(keep);

% ---- unique regions in the requested order ----
regionList = unique(string(patMap(:,2)),'stable');

for r = 1:numel(regionList)
    reg = regionList(r);
    T   = combined_region(strcmp(combined_region.Region, reg), :);
    if isempty(T)
        warning('No rows found for region "%s". Skipping.', reg);
        continue;
    end

    % --- collapse across days for this region (same logic as single-plot version) ---
    mu    = nan(M,1);
    CI    = nan(2,M);
    pComb = nan(M,1);
    nPer  = nan(M,1);

    for i = 1:M
        sel   = strcmp(T.Pair, pairOrder{i});
        vals  = T.(effCol)(sel);
        pvals = T.(pCol)(sel);

        vals  = vals(isfinite(vals));
        pvals = pvals(isfinite(pvals));

        nPer(i) = numel(vals);
        if isempty(vals), continue; end

        mu(i) = median(vals,'omitnan');

        if numel(vals)==1
            CI(:,i) = [vals; vals];
        else
            bs = nan(NBoot,1);
            for b = 1:NBoot
                ii = randi(numel(vals), numel(vals), 1);
                bs(b) = median(vals(ii),'omitnan');
            end
            CI(:,i) = prctile(bs,[2.5 97.5])';
        end

        if ~isempty(pvals)
            pvals = max(pvals, 1e-12);    % floor to avoid -Inf
            X  = -2*sum(log(pvals));
            df = 2*numel(pvals);
            pComb(i) = 1 - chi2cdf(X, df);
        end
    end

    % --- plotting for this region ---
    figure('Color','w','Name',sprintf('Region %s (%s)', reg, char(metric)));
    hold on;
    x = (1:M).';

    if showPts
        gcol = 0.7*[1 1 1];
        for i = 1:M
            sel  = strcmp(T.Pair, pairOrder{i});
            vals = T.(effCol)(sel);
            vals = vals(isfinite(vals));
            if isempty(vals), continue; end
            jitter = (rand(size(vals))-0.5)*0.15;
            plot(i + jitter, vals, '.', 'Color', gcol, 'MarkerSize', 24);
        end
    end

    mu_col  = mu(:);
    ci_low  = CI(1,:).';
    ci_high = CI(2,:).';
    yneg    = mu_col - ci_low;
    ypos    = ci_high - mu_col;

    errorbar(x, mu_col, yneg, ypos, 'o', 'LineWidth',1.8, 'CapSize',10, 'MarkerSize',6);

    finiteCI = isfinite(ci_low) & isfinite(ci_high);
    rngY = diff([min(ci_low(finiteCI)), max(ci_high(finiteCI))]);
    if isempty(rngY) || ~isfinite(rngY) || rngY==0, rngY = 0.1; end

    for i = 1:M
        if ~isfinite(pComb(i)), continue; end
        if     pComb(i) < 0.001, stars='***';
        elseif pComb(i) < 0.01,  stars='**';
        elseif pComb(i) < alpha, stars='*';
        else,  stars='n.s.';
        end
        yTop = ci_high(i); if ~isfinite(yTop), yTop = mu_col(i); end
        text(x(i), yTop + 0.03*rngY, stars, 'HorizontalAlignment','center','FontWeight','bold');
    end

    xlim([0.5 M+0.5]); ylim([-1 1])
    xticks(x); xticklabels(pairDisplay);
    ylabel(ylab);
    title(sprintf('Region: %s — median with bootstrap 95%% CI; p via Fisher''s method', reg));
    grid on; box off;

    annotation('textbox',[0.13 0.01 0.8 0.05], 'String', ...
        sprintf('Points = per-day estimates. Combined p via Fisher. n per pair: %s', ...
        strjoin(arrayfun(@(n)sprintf('%d',n), nPer,'UniformOutput',false), ', ')), ...
        'EdgeColor','none','HorizontalAlignment','left','FontSize',14);
end
end

% ---------- helper ----------
function s = pformat(p)
    if ~isfinite(p) || isnan(p), s = 'NaN'; return; end
    if p < 1e-3
        s = sprintf('%.2e', p);
    else
        s = sprintf('%.3f', p);
    end
end
