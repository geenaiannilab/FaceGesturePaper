%% === Gather chi-square results into a single table ===
% Assumes variables from your snippet already exist in workspace:
% pThres, alpha, regions, cortRegionOut, and all tbl_*, chi2stat_*, pval_* variables.

% Helper to compute df from crosstab table
alpha = 0.05
getDF = @(T) (size(T,1)-1) * (size(T,2)-1);

% Container for rows
Stats = table();
ContingencyTables = {};

addRow = @(cmp,factor,regionsCompared,chi2,p,alpha,df,tbl) struct( ...
    'Comparison',        string(cmp), ...
    'Factor',            string(factor), ...
    'RegionsCompared',   string(regionsCompared), ...
    'ChiSquare',         chi2, ...
    'DF',                df, ...
    'PValue',            p, ...
    'AlphaUsed',         alpha, ...
    'Significant',       p < alpha, ...
    'Table',             {tbl} );

% --- Proper empty struct initialization (prevents "conversion to double" error)
rows = struct( ...
    'Comparison',      string.empty(0,1), ...
    'Factor',          string.empty(0,1), ...
    'RegionsCompared', string.empty(0,1), ...
    'ChiSquare',       double.empty(0,1), ...
    'DF',              double.empty(0,1), ...
    'PValue',          double.empty(0,1), ...
    'AlphaUsed',       double.empty(0,1), ...
    'Significant',     false(0,1), ...
    'Table',           cell(0,1) );


%% ---- Overall tests across all 4 regions ----
rows(end+1) = addRow('Overall', 'Behavior', 'S1,M1,PMv,M3', chi2statBhv, pvalBhv, pThres, getDF(tblBhv), tblBhv);
rows(end+1) = addRow('Overall', 'Time',     'S1,M1,PMv,M3', chi2statTime, pvalTime, pThres, getDF(tblTime), tblTime);
rows(end+1) = addRow('Overall', 'Interact', 'S1,M1,PMv,M3', chi2statInt,  pvalInt,  pThres, getDF(tblInt),  tblInt);

%% ---- Pairwise TIME post-hoc (alpha = alpha) ----
rows(end+1) = addRow('Pairwise','Time','S1 vs M1', chi2stat_S1M1time, pval_S1M1time, alpha, getDF(tbl_S1M1time), tbl_S1M1time);
rows(end+1) = addRow('Pairwise','Time','S1 vs PMv', chi2stat_S1PMvtime, pval_S1PMvtime, alpha, getDF(tbl_S1PMvtime), tbl_S1PMvtime);
rows(end+1) = addRow('Pairwise','Time','S1 vs M3', chi2stat_S1M3time, pval_S1M3time, alpha, getDF(tbl_S1M3time), tbl_S1M3time);
rows(end+1) = addRow('Pairwise','Time','M1 vs PMv', chi2stat_M1PMvtime, pval_M1PMvtime, alpha, getDF(tbl_M1PMvtime), tbl_M1PMvtime);
rows(end+1) = addRow('Pairwise','Time','M1 vs M3', chi2stat_M1M3time, pval_M1M3time, alpha, getDF(tbl_M1M3time), tbl_M1M3time);
rows(end+1) = addRow('Pairwise','Time','PMv vs M3', chi2stat_PMvM3time, pval_PMvM3time, alpha, getDF(tbl_PMvM3time), tbl_PMvM3time);

%% ---- Pairwise INTERACTION post-hoc (alpha = alpha) ----
rows(end+1) = addRow('Pairwise','Interact','S1 vs M1', chi2stat_S1M1int, pval_S1M1int, alpha, getDF(tbl_S1M1int), tbl_S1M1int);
rows(end+1) = addRow('Pairwise','Interact','S1 vs PMv', chi2stat_S1PMvint, pval_S1PMvint, alpha, getDF(tbl_S1PMvint), tbl_S1PMvint);
rows(end+1) = addRow('Pairwise','Interact','S1 vs M3', chi2stat_S1M3int, pval_S1M3int, alpha, getDF(tbl_S1M3int), tbl_S1M3int);
rows(end+1) = addRow('Pairwise','Interact','M1 vs PMv', chi2stat_M1PMvint, pval_M1PMvint, alpha, getDF(tbl_M1PMvint), tbl_M1PMvint);
rows(end+1) = addRow('Pairwise','Interact','M1 vs M3', chi2stat_M1M3int, pval_M1M3int, alpha, getDF(tbl_M1M3int), tbl_M1M3int);
rows(end+1) = addRow('Pairwise','Interact','PMv vs M3', chi2stat_PMvM3int, pval_PMvM3int, alpha, getDF(tbl_PMvM3int), tbl_PMvM3int);

%% Convert to table (stats) + cell array (contingency tables)
Stats = struct2table(rows);
ContingencyTables = Stats.Table;           % cell array, same row order
Stats.Table = [];                          % drop from the stats table for CSV export

%% (Optional) Add adjusted p-value column for clarity
Stats.PValueAdj = nan(height(Stats),1);
isPair = Stats.Comparison=="Pairwise";
Stats.PValueAdj(isPair) = Stats.PValue(isPair);  % equals raw p but judged vs alpha
Stats = movevars(Stats, 'PValueAdj', 'After', 'PValue');

%% Display quick preview
disp(Stats);

%% Save outputs
writetable(Stats, 'chi2_summary_stats.csv');        % portable summary (no matrices)
save('/Users/geena/Dropbox/PhD/SUAinfo/chi2_summary_full.mat','Stats','ContingencyTables');  % full, includes crosstab tables

%% === Append per-region proportions to the Overall rows ===
% Prereqs: 'regions' is a cellstr or string array with region names in the same
% order as the coded IDs in 'cortRegionOut' (1..numel(regions)).
% Uses crosstab outputs: tblBhv, tblTime, tblInt (rows = [0;1], cols = regions)

assert(exist('regions','var')==1, 'Variable "regions" not found.');
R = numel(regions);

%% ---- Robust region names aligned to crosstab column order ----
regIDs = sort(unique(cortRegionOut(:)'));   % e.g., [1 2 3 4]

% Default names (safe fallback)
if isequal(regIDs, [1 2 3 4])
    defaultNames = ["S1","M1","PMv","M3"]';
else
    defaultNames = arrayfun(@(k) "Region"+k, regIDs, 'UniformOutput', true).';
end
regNamesRaw = defaultNames;

% Try to use user-provided "regions" if present, supporting many shapes
if exist('regions','var')
    cand = strings(0,1);
    try
        if isstring(regions)
            cand = regions(:);

        elseif iscellstr(regions)          % cell of char
            cand = string(regions(:));

        elseif iscell(regions)             % heterogeneous cell
            cand = strings(numel(regions),1);
            for i = 1:numel(regions)
                x = regions{i};
                if ischar(x) || isstring(x)
                    cand(i) = string(x);
                elseif isnumeric(x) && isscalar(x)
                    cand(i) = "Region" + string(x);
                elseif isstruct(x)
                    if isfield(x,'name'),   cand(i) = string(x.name);
                    elseif isfield(x,'label'), cand(i) = string(x.label);
                    else, cand(i) = missing;
                    end
                else
                    cand(i) = missing;
                end
            end

        elseif isstruct(regions)           % struct array with name/label
            if isfield(regions,'name')
                cand = string({regions.name}.');
            elseif isfield(regions,'label')
                cand = string({regions.label}.');
            end

        elseif istable(regions)            % table with Name/Label column
            vn = regions.Properties.VariableNames;
            if any(strcmpi(vn,'Name'))
                cand = string(regions{:,'Name'});
            elseif any(strcmpi(vn,'Label'))
                cand = string(regions{:,'Label'});
            end
        end
    catch
        % swallow and keep fallback
    end

    % Only accept if the size matches and entries are nonempty
    if numel(cand) == numel(regIDs) && all(strlength(cand) > 0)
        regNamesRaw = cand(:);
    end
end

% Final: valid MATLAB identifiers for table column headers
regNames = matlab.lang.makeValidName(regNamesRaw, 'ReplacementStyle','delete');


% Identify the three Overall rows (one each for Behavior, Time, Interact)
maskOverall = Stats.Comparison == "Overall";

% Behavior overall row gets per-region stats from tblBhv
maskBhv = maskOverall & Stats.Factor == "Behavior";
if any(maskBhv)
    Stats = add_region_cols(Stats, maskBhv, tblBhv, regNames);
end

% Time overall row from tblTime
maskTime = maskOverall & Stats.Factor == "Time";
if any(maskTime)
    Stats = add_region_cols(Stats, maskTime, tblTime, regNames);
end

% Interaction overall row from tblInt
maskInt = maskOverall & Stats.Factor == "Interact";
if any(maskInt)
    Stats = add_region_cols(Stats, maskInt, tblInt, regNames);
end

%% Re-save with the new columns
writetable(Stats, 'chi2_summary_stats.csv');
save('/Users/geena/Dropbox/PhD/SUAinfo/chi2_summary_full.mat','Stats','ContingencyTables');

%% Preview the Overall rows with the new columns
disp(Stats(maskOverall,:));

%% === Add Benjamini–Hochberg FDR (per factor) to pairwise rows ===
% Prereqs: Stats table exists with columns: Comparison, Factor, PValue.
% Uses alpha = pThres for the BH level (change if desired).

% Initialize columns
if ~ismember('q_BH', string(Stats.Properties.VariableNames))
    Stats.q_BH   = nan(height(Stats),1);
end
if ~ismember('Sig_BH', string(Stats.Properties.VariableNames))
    Stats.Sig_BH = false(height(Stats),1);
end
if ~ismember('BH_CritP', string(Stats.Properties.VariableNames))
    Stats.BH_CritP = nan(height(Stats),1);
end

pairMask = (Stats.Comparison == "Pairwise");
alphaBH  = pThres;   % set your desired FDR level here (e.g., 0.05)

for fac = ["Time","Interact"]
    m = find(pairMask & Stats.Factor == fac);
    if isempty(m), continue; end
    pvec = Stats.PValue(m);

    % Run BH FDR (helper defined below or put in its own file)
    [qvals, sigMask, crit_p] = bh_fdr_vector(pvec, alphaBH);

    % Write back
    Stats.q_BH(m)     = qvals;
    Stats.Sig_BH(m)   = sigMask;
    Stats.BH_CritP(m) = crit_p;   % same value repeated for the family
end

% Re-save with new columns
writetable(Stats, 'chi2_summary_stats.csv');
save('chi2_summary_full.mat','Stats','ContingencyTables');

%% === Print APA-style lines for each pairwise chi-square test (Unicode symbols) ===
pairMask = (Stats.Comparison == "Pairwise");
reportLines = strings(sum(pairMask),1);
ii = find(pairMask);

for k = 1:numel(ii)
    r = ii(k);
    tbl     = ContingencyTables{r};
    N       = sum(tbl(:));
    factor  = char(Stats.Factor(r));            % "Time" or "Interact"
    compare = char(Stats.RegionsCompared(r));   % e.g., "S1 vs M1"
    df      = Stats.DF(r);
    chi2val = Stats.ChiSquare(r);
    pval    = Stats.PValue(r);
    alphaUse= Stats.AlphaUsed(r);

    % APA-style p text (no leading zero; < .001 when appropriate)
    if pval < 0.001
        pText = 'p < .001';
    else
        pText = sprintf('p = %.3f', pval);
        pText = regexprep(pText,'p = 0\.','p = .');  % remove leading zero only
    end

    % Significance text
    if Stats.Significant(r)
        sigTxt = 'significant';
    else
        sigTxt = 'ns';
    end

    % Build Unicode string
    line = sprintf('%s, %s: χ²(%d, N=%d) = %.2f, %s, %s (α_adj = %.4g)', ...
                   factor, compare, df, N, chi2val, pText, sigTxt, alphaUse);

    reportLines(k) = string(line);
end

% Print to Command Window
disp("Pairwise chi-square tests (APA-style, Unicode):");
for k = 1:numel(reportLines)
    disp(reportLines(k));
end


% Also write to a text file
fid = fopen('chi2_pairwise_report.txt','w','n','UTF-8'); % ensure UTF-8 encoding
fprintf(fid, 'Pairwise chi-square tests (APA-style, Unicode):\n');
for k = 1:numel(reportLines)
    fprintf(fid, '%s\n', reportLines(k));
end
fclose(fid);



function [qvals, sigMask, crit_p] = bh_fdr_vector(p_in, alpha)
% bh_fdr_vector  Benjamini–Hochberg FDR for a numeric vector (NaNs allowed).
% INPUTS: p_in (nx1 or 1xn), alpha (e.g., 0.05).  OUTPUTS: qvals (size p_in),
% sigMask (logical, same size), crit_p (scalar BH critical p or NaN).

    if nargin < 2 || isempty(alpha), alpha = 0.05; end
    p_in = p_in(:);
    valid = isfinite(p_in);
    p = p_in(valid);

    qvals = nan(size(p_in));
    sigMask = false(size(p_in));
    crit_p = NaN;

    if isempty(p), return; end

    [ps, order] = sort(p,'ascend');
    m = numel(ps);
    ranks = (1:m)';

    % q-values (BH)
    q = (m ./ ranks) .* ps;
    q = min(1, cummin(q(end:-1:1)));
    q = q(end:-1:1);

    % Discoveries
    thresh = (ranks / m) * alpha;
    pass = ps <= thresh;
    if any(pass)
        k = find(pass, 1, 'last');
        crit_p = ps(k);
        sig_sorted = false(m,1);
        sig_sorted(1:k) = true;
    else
        sig_sorted = false(m,1);
        crit_p = NaN;
    end

    % Map back
    q_back = nan(m,1);         q_back(order) = q;
    sig_back = false(m,1);     sig_back(order) = sig_sorted;

    qvals(valid) = q_back;
    sigMask(valid) = sig_back;
end

% Function: add N_<Region> and PropSig_<Region> columns from an overall crosstab
function Stats = add_region_cols(Stats, rowMask, tbl, regNames)
    % Detect TRUE row (significant=1). crosstab on logical typically [0;1]
    % Fall back to "take last row" if size > 1.
    if size(tbl,1) >= 2
        trueRow = size(tbl,1); % robust even if categories weren’t logical
    else
        error('Crosstab table does not have at least two rows (expected [notSig; sig]).');
    end

    % Column-wise totals and sig counts
    N_by_reg   = sum(tbl,1);           % total per region
    Sig_by_reg = tbl(trueRow, :);      % significant per region
    Prop_by_reg = Sig_by_reg ./ max(N_by_reg, 1);  % avoid divide-by-zero

    % Ensure columns exist (add if missing), then fill only the masked rows
    for c = 1:numel(regNames)
        nCol  = "N_" + regNames(c);
        pCol  = "PropSig_" + regNames(c);
        if ~ismember(nCol, string(Stats.Properties.VariableNames))
            Stats.(nCol) = NaN(height(Stats),1);
        end
        if ~ismember(pCol, string(Stats.Properties.VariableNames))
            Stats.(pCol) = NaN(height(Stats),1);
        end
        Stats.(nCol)(rowMask) = N_by_reg(c);
        Stats.(pCol)(rowMask) = Prop_by_reg(c);
    end
end
