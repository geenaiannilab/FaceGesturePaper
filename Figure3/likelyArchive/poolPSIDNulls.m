%%%% written 250201 
% Pool kinematicDecoding and kinematicDecodingStd across 25 iterations where:
%   - Base file contains iterations 1 & 2 as columns 1 and 2
%   - Each subsequent file N=j contains data in column j (cols < j are zeros)
%
% Output:
%   agg.combined.(region).kinematicDecoding     -> [nRows x 25]
%   agg.combined.(region).kinematicDecodingStd  -> [nRows x 25]
%   agg.info.* with bookkeeping
%
% Edit CONFIG to your paths/names if needed.

% ----------------------- CONFIG -----------------------
pseudoPop  = 'LSonly/balanced';
subj       = 'combinedSubjs';
popSize    = '50';
nIters     = 25;

rootDir    = '/Users/geena/Dropbox/PhD/SUAinfo/PSIDResults2023/Pseudopopulations';
resultsDir = fullfile(rootDir, pseudoPop, subj, ['N_' popSize 'cells'], 'results');

baseFile   = fullfile(resultsDir, 'copy_PP1-2_decodingNullResults.mat');     % iter 1 & 2
filePat    = 'decodingNullResults_pp_N%d.mat';                           % iter 3..25
outFile    = fullfile(resultsDir, 'decodingNullResults_pp_pooled_byCol.mat');

fieldsToPool = {'kinematicDecoding','kinematicDecodingStd'};
% ------------------------------------------------------

% ----- load base to get regions and seed columns 1 & 2 -----
Sbase = load(baseFile);
assert(isfield(Sbase,'shuffledResultsAll') && isfield(Sbase.shuffledResultsAll,'combined'), ...
    'Base file must contain shuffledResultsAll.combined.*');
Cbase   = Sbase.shuffledResultsAll.combined;
regions = fieldnames(Cbase);

% init containers
agg = struct('combined', struct());
for r = 1:numel(regions)
    reg = regions{r};
    for f = 1:numel(fieldsToPool)
        fld = fieldsToPool{f};
        agg.combined.(reg).(fld) = []; % will become [nRows x nIters]
    end
end

% ----- seed from base file: take col 1 and col 2 if present -----
for r = 1:numel(regions)
    reg = regions{r};
    for f = 1:numel(fieldsToPool)
        fld = fieldsToPool{f};
        if ~isfield(Cbase.(reg), fld) || isempty(Cbase.(reg).(fld)), continue; end

        M = Cbase.(reg).(fld);
        c1 = extract_col(M, 1, baseFile, reg, fld);   % iteration 1
        c2 = extract_col(M, 2, baseFile, reg, fld);   % iteration 2

        % Create [nRows x nIters] and drop in the first two columns
        nRows = max([numel(c1), numel(c2), 0]);
        if nRows == 0
            agg.combined.(reg).(fld) = nan(0, nIters);
            continue
        end
        A = nan(nRows, nIters, 'like', double(1));

        if ~isempty(c1), A(:,1) = pad_or_trim(c1(:), nRows); end
        if ~isempty(c2), A(:,2) = pad_or_trim(c2(:), nRows); end
        agg.combined.(reg).(fld) = A;
    end
end

% ----- loop over iteration files 3..nIters and take column j -----
for j = 3:nIters
    fn = fullfile(resultsDir, sprintf(filePat, j));
    if ~isfile(fn)
        warning('Missing file: %s (skipping iter %d)', fn, j);
        continue
    end
    S = load(fn);
    if ~isfield(S,'shuffledResultsAll') || ~isfield(S.shuffledResultsAll,'combined')
        warning('No shuffledResultsAll.combined in %s (skipping iter %d)', fn, j);
        continue
    end
    C = S.shuffledResultsAll.combined;

    for r = 1:numel(regions)
        reg = regions{r};
        if ~isfield(C, reg)
            warning('Region %s missing in %s (skipping)', reg, fn);
            continue
        end
        for f = 1:numel(fieldsToPool)
            fld = fieldsToPool{f};
            if ~isfield(C.(reg), fld) || isempty(C.(reg).(fld)), continue; end

            colj = extract_col(C.(reg).(fld), j, fn, reg, fld);   % column j only
            if isempty(colj), continue; end

            % assign into agg.combined.(reg).(fld)(:, j), padding/trimming rows
            A = agg.combined.(reg).(fld);
            nRows = size(A,1);
            cj    = pad_or_trim(colj(:), nRows);
            if numel(cj) ~= nRows
                % First time encountering this field with different row count (e.g., if base had 0 rows)
                nRowsNew = numel(cj);
                if nRows == 0
                    A = nan(nRowsNew, nIters, 'like', cj);
                else
                    % Expand existing A to the new row count
                    B = nan(nRowsNew, nIters, 'like', A);
                    B(1:nRows, :) = A;
                    A = B;
                end
                nRows = size(A,1);
                cj    = pad_or_trim(colj(:), nRows);
            end
            A(:, j) = cj;
            agg.combined.(reg).(fld) = A;
        end
    end
end

% ----- annotate + save -----
agg.info.resultsDir = resultsDir;
agg.info.baseFile   = baseFile;
agg.info.nIters     = nIters;
save(outFile, 'agg', '-v7.3');
fprintf('Saved pooled results to:\n  %s\n', outFile);


% ====================== helpers ======================

function col = extract_col(M, j, fn, reg, fld)
% Extract exactly column j (iteration j) from M, tolerating trailing dims.
% If M has fewer than j columns in dim-2, return [] with a warning.
if ndims(M) == 2
    if size(M,2) < j
        warning('%s.%s in %s has only %d columns; need col %d.', reg, fld, fn, size(M,2), j);
        col = [];
        return
    end
    col = M(:, j);
elseif ndims(M) >= 3
    if size(M,2) < j
        warning('%s.%s in %s has only %d columns; need col %d.', reg, fld, fn, size(M,2), j);
        col = [];
        return
    end
    tmp = squeeze(M(:, j, :));   % collapse trailing dims
    if isvector(tmp)
        col = tmp(:);
    elseif size(tmp,2) == 1
        col = tmp;               % [nRows x 1]
    else
        % Ambiguous trailing dims; take first trailing column to proceed
        warning('%s.%s in %s has trailing dims; using first trailing column for col %d.', reg, fld, fn, j);
        col = tmp(:,1);
    end
else
    col = [];
end
end

function v = pad_or_trim(v, T)
% Force v to length T by padding with NaN or truncating.
v = v(:);
n = numel(v);
if n < T
    v = [v; nan(T-n,1)];
elseif n > T
    v = v(1:T);
end
end
