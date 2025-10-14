
close all; clear all;

load('Matfiles/markerTrajectoriesForMovies.mat');

nBootstrap = 100; 

nThreats = size(bhvTraj(1).data,2);
nLS = size(bhvTraj(2).data,2);
nChew = size(bhvTraj(3).data,2);
nTrials = [nThreats nLS nChew];
colorMap = [1 0 0 ; 0 0 1; 0 1 0];
Time = taxis*1000;
msPerFrame = 10;
onsetIdx = find(taxis == 0);

outDir = '~/Desktop/';

% Call the function
make_marker_movies_biggestTrial(markerData2plot, ...
    'Time', Time, ...
    'MsPerFrame', msPerFrame, ...
    'OnsetIdx', onsetIdx, ...
    'NTrail', 5, ...
    'Colors', [1 0 0; 0 0 1; 0 1 0], ...
    'GestureNames', {'threat','lipsmack','chew'}, ...
    'InvertY', true, ...
    'FPS', 24, ...
    'FrameStep', 1, ...
    'OutDir', outDir);

% Compute (diag Σ is usually fine in PC space)
out = mahal_from_trials(bhvTraj, 'CovMode','diag', 'B', nBootstrap, 'Reg', 1e-6, 'Verbose', true);

% Plot: thick mean + shaded ±2·SEM
plot_mahal_mean_sem(out, taxis);

% quick plot
t = (1:size(bhvTraj(1).avg,1))';
figure; 
plot(taxis, out.D,'linew',8); legend({'Thr-LS', 'Thr-Ch', 'LS-Ch'}, 'Location','best'); xlabel('Time'); ylabel('Mahalanobis distance');
title('Pairwise Mahalanobis distances between gesture trajectories');

%%
make_marker_movies(markerData2plot, ...
    'Time', taxis, ...
    'TopFrac', 0.01, 'MinTrials', 5, ...
    'OutDir', './movies', 'FPS', 48, 'FrameStep', 1, ...
    'InvertY', true, 'Colors', colorMap);

make_marker_movies_meanfade(markerData2plot, ...
    'Time', taxis, ...             % [T x 1] with 0 = onset
    'TopFrac', 0.05, 'MinTrials', 3, ...
    'NTrail', 15, ...                 % length of fading tail
    'FPS',48,...
    'FrameStep',2,...
    'Colors', colorMap, ...
    'GestureNames', {'threat','lipsmack','chew'}, ...
    'OutDir','./');

function out = mahal_trajectories(bhvTraj, varargin)

% MAHAL_TRAJECTORIES
% Compute time-resolved Mahalanobis distances between three gesture
% trajectories in PC space from bhvTraj(1x3) with fields:
%   .avg [T x P]  : mean PC trajectory
%   .sem [T x P]  : SEM of PC trajectory
% Optional:
%   .cov [P x P x T] : (per gesture) full covariance per time (if available)
%
% out.D   : [T x 3] distances for pairs {'1-2','1-3','2-3'}
% out.pairs : pair labels
% out.distMatPerTime : 1xT cell, each a 3x3 matrix of pairwise distances
% out.opts : resolved options

p = inputParser;
addParameter(p, 'CovMode', 'pooled-diag', @(s)ischar(s)||isstring(s)); % 'pooled-diag' | 'pooled-full' | 'from-cov-field'
addParameter(p, 'Trials', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==3)); % [n1 n2 n3]
addParameter(p, 'Reg', 1e-6, @(x)isnumeric(x)&&isscalar(x)&&x>=0);  % ridge added to diag/Σ
addParameter(p, 'TimeIdx', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
parse(p, varargin{:});
opts = p.Results;

% Basic checks
assert(numel(bhvTraj)==3, 'Provide bhvTraj as a 1x3 struct (one per gesture).');
T = size(bhvTraj(1).avg,1);
P = size(bhvTraj(1).avg,2);
for k = 1:3
    assert(isequal(size(bhvTraj(k).avg), [T P]), 'avg size mismatch for gesture %d', k);
    assert(isequal(size(bhvTraj(k).sem), [T P]), 'sem size mismatch for gesture %d', k);
end
if isempty(opts.TimeIdx), opts.TimeIdx = 1:T; end
TI = opts.TimeIdx(:)';

pairs = [1 2; 1 3; 2 3];
pairLabels = {'1-2','1-3','2-3'};

D = nan(T,3);
distMatPerTime = cell(1,T);

% Helper to form pooled Σ at a time t for a pair (i,j)
    function Sigma = makeSigma(t,i,j)
        switch lower(opts.CovMode)
            case 'pooled-diag'
                % Back out variances from SEM if Trials provided
                if ~isempty(opts.Trials)
                    n_i = opts.Trials(i); n_j = opts.Trials(j);
                    var_i = (bhvTraj(i).sem(t,:)).^2 .* max(n_i,1);
                    var_j = (bhvTraj(j).sem(t,:)).^2 .* max(n_j,1);
                    % classic pooled variance per dimension
                    var_p = ((max(n_i-1,1)).*var_i + (max(n_j-1,1)).*var_j) ./ max(n_i+n_j-2,1);
                else
                    % No trial counts: use average of SEM^2 as a proxy scale
                    var_p = 0.5*((bhvTraj(i).sem(t,:)).^2 + (bhvTraj(j).sem(t,:)).^2);
                end
                % Regularize
                var_p = var_p + opts.Reg;
                Sigma = diag(var_p);

            case 'pooled-full'
                % Estimate full Σ from SEMs assuming diagonal (fallback) then inflate by a correlation guess of zero.
                % If you truly need full Σ, provide bhvTraj(k).cov(:,:,t) and use 'from-cov-field'.
                warning_once('Using ''pooled-full'' without cov fields; falling back to diagonal.');
                var_p = 0.5*((bhvTraj(i).sem(t,:)).^2 + (bhvTraj(j).sem(t,:)).^2) + opts.Reg;
                Sigma = diag(var_p);

            case 'from-cov-field'
                assert(isfield(bhvTraj(i),'cov') && ~isempty(bhvTraj(i).cov), ...
                       'from-cov-field requires bhvTraj(%d).cov(:,:,t)', i);
                assert(isfield(bhvTraj(j),'cov') && ~isempty(bhvTraj(j).cov), ...
                       'from-cov-field requires bhvTraj(%d).cov(:,:,t)', j);
                Ci = bhvTraj(i).cov(:,:,t);
                Cj = bhvTraj(j).cov(:,:,t);
                assert(isequal(size(Ci), [P P]) && isequal(size(Cj), [P P]), 'cov size mismatch at t=%d', t);
                % Pooled full covariance
                if ~isempty(opts.Trials)
                    n_i = opts.Trials(i); n_j = opts.Trials(j);
                    Sigma = ((max(n_i-1,1))*Ci + (max(n_j-1,1))*Cj) / max(n_i+n_j-2,1);
                else
                    Sigma = 0.5*(Ci + Cj);
                end
                % Regularize
                Sigma = Sigma + opts.Reg*eye(P);
            otherwise
                error('Unknown CovMode: %s', opts.CovMode);
        end
    end

% Optional one-time warning control
function warning_once(msg)
    persistent didWarn
    if isempty(didWarn) || ~didWarn
        warning(msg);
        didWarn = true;
    end
end

% Main loop over time
for t = TI
    dm = zeros(3,3); % 3x3 matrix of distances at this time
    for pidx = 1:3
        i = pairs(pidx,1); j = pairs(pidx,2);
        mu_i = bhvTraj(i).avg(t,:); 
        mu_j = bhvTraj(j).avg(t,:);
        delta = (mu_i - mu_j).';
        Sigma = makeSigma(t,i,j);

        % Mahalanobis distance
        % Use Cholesky solve if possible for stability
        [L,pchol] = chol(Sigma,'lower');
        if pchol==0
            y = L \ delta;
            d2 = y.'*y;  % (L^{-1} delta)'(L^{-1} delta) == delta' Sigma^{-1} delta
        else
            % Fallback to pinv if Sigma not PD
            iS = pinv(Sigma);
            d2 = delta.'*iS*delta;
        end
        d = sqrt(max(d2,0));
        D(t,pidx) = d;
        dm(i,j) = d; dm(j,i) = d;
    end
    distMatPerTime{t} = dm; % diagonal stays 0
end

out = struct();
out.D = D;                     % [T x 3], columns correspond to pairLabels
out.pairs = pairLabels;
out.distMatPerTime = distMatPerTime;
out.opts = opts;
end

function out = mahal_from_trials(bhvTraj, varargin)

% Compute pairwise Mahalanobis distance time-courses between 3 gestures
% from per-trial PC trajectories, with bootstrap SEM bands.
%
% INPUT (per gesture g = 1..3)
%   bhvTraj(g).data : [T x N_g x P]  (time x trials x PCs)
%
% OPTIONS
%   'CovMode' : 'diag' | 'full'  (default 'diag'; 'full' needs enough trials)
%   'B'       : bootstrap resamples per gesture (default 500)
%   'Reg'     : ridge regularization added to Σ (default 1e-6)
%   'Pairs'   : which pairs to compute (default [1 2; 1 3; 2 3])
%   'Verbose' : true/false (default false)
%
% OUTPUT
%   out.t             : [T x 1] time index (1..T)
%   out.pairs         : labels
%   out.D             : [T x 3] point-estimate distances from all trials
%   out.Dsem          : [T x 3] bootstrap SEM across resamples
%   out.Dboot         : [T x 3 x B] bootstrap sample distances
%   out.opts          : resolved options

p = inputParser;
addParameter(p,'CovMode','diag',@(s)ischar(s)||isstring(s));
addParameter(p,'B',100,@(x)isnumeric(x)&&isscalar(x)&&x>=10);
addParameter(p,'Reg',1e-6,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'Pairs',[1 2; 1 3; 2 3],@(x)isnumeric(x)&&size(x,2)==2&&all(x(:)>=1&x(:)<=3));
addParameter(p,'Verbose',false,@islogical);
parse(p,varargin{:});
opt = p.Results;

assert(numel(bhvTraj)==3, 'Expect 1x3 structure array.');

% Basic sizes
[T, N1, P] = size(bhvTraj(1).data); %#ok<ASGLU>
for g = 1:3
    sz = size(bhvTraj(g).data);
    assert(sz(1)==T && sz(3)==P, 'All gestures must share T and P.');
end

pairs = opt.Pairs;
pairLabels = arrayfun(@(k)sprintf('%d-%d',pairs(k,1),pairs(k,2)), 1:size(pairs,1), 'uni',0);
K = size(pairs,1);

% ---------- helpers ----------
    function [mu, C, n] = mean_cov_at_t(g, t, mode, reg)
        % X: [N x P] trials-by-PCs at time t for gesture g
        Xi = squeeze(bhvTraj(g).data(t,:,:));   % [N x P]
        if isvector(Xi), Xi = Xi(:)'; end
        good = all(isfinite(Xi),2);
        Xi = Xi(good,:); n = size(Xi,1);
        mu = mean(Xi, 1);
        if strcmpi(mode,'full')
            if n>=2
                C = cov(Xi, 1);   % MLE (divide by N), better conditioned
            else
                C = eye(P).*eps;
            end
        else % 'diag'
            if n>=1
                v = var(Xi, 1, 1); % [1 x P]
                C = diag(v);
            else
                C = eye(P).*eps;
            end
        end
        C = C + reg*eye(P);
    end

    function d = maha_pair_at_t(i,j,t)
        [mui, Ci, ni] = mean_cov_at_t(i,t,opt.CovMode,opt.Reg);
        [muj, Cj, nj] = mean_cov_at_t(j,t,opt.CovMode,opt.Reg);
        delta = (mui - muj).';
        % pooled Σ
        if ni+nj>=3
            Sig = ((max(ni-1,1))*Ci + (max(nj-1,1))*Cj) / max(ni+nj-2,1);
        else
            Sig = 0.5*(Ci + Cj);
        end
        % Mahalanobis via Cholesky
        [L,pchol] = chol(Sig,'lower');
        if pchol==0
            y = L \ delta;
            d2 = y.'*y;
        else
            d2 = max(delta.'*(pinv(Sig))*delta, 0);
        end
        d = sqrt(d2);
    end
% -----------------------------

% Point estimate from all trials
D = nan(T,K);
for t = 1:T
    for k = 1:K
        D(t,k) = maha_pair_at_t(pairs(k,1), pairs(k,2), t);
    end
end

% Bootstrap SEM (resample trials WITHIN gesture independently)
B = opt.B;
Dboot = nan(T,K,B);
Ns = arrayfun(@(g) size(bhvTraj(g).data,2), 1:3);
if opt.Verbose, fprintf('Bootstrapping %d resamples...', B); end
for b = 1:B
    % draw indices
    idx = cell(1,3);
    for g = 1:3
        n = Ns(g);
        idx{g} = randi(n, [n,1]);  % [n x 1], with replacement
    end
    % temporary struct with resampled trials
    bhvRes = bhvTraj;
    for g = 1:3
        bhvRes(g).data = bhvTraj(g).data(:, idx{g}, :);
    end

    % compute distances for this bootstrap sample
    for t = 1:T
        for k = 1:K
            i = pairs(k,1); j = pairs(k,2);

            [mui,Ci,ni] = local_mcov(bhvRes,i,t,opt.CovMode,opt.Reg);
            [muj,Cj,nj] = local_mcov(bhvRes,j,t,opt.CovMode,opt.Reg);
            delta = (mui - muj).';
            if ni+nj>=3
                Sig = ((max(ni-1,1))*Ci + (max(nj-1,1))*Cj) / max(ni+nj-2,1);
            else
                Sig = 0.5*(Ci + Cj);
            end
            [L,pchol] = chol(Sig,'lower');
            if pchol==0
                y = L \ delta;
                d2 = y.'*y;
            else
                d2 = max(delta.'*(pinv(Sig))*delta, 0);
            end
            Dboot(t,k,b) = sqrt(d2);
        end
    end
end
if opt.Verbose, fprintf(' done.\n'); end

Dsem = std(Dboot, 0, 3);   % bootstrap SE

out = struct();
out.t = (1:T).';
out.pairs = pairLabels;
out.D = D;
out.Dsem = Dsem;
out.Dboot = Dboot;
out.opts = opt;

    function [mu,C,n] = local_mcov(S,g,t,mode,reg)
        X = squeeze(S(g).data(t,:,:)); if isvector(X), X = X(:)'; end
        good = all(isfinite(X),2); X = X(good,:); n = size(X,1);
        mu = mean(X,1);
        if strcmpi(mode,'full')
            if n>=2, C = cov(X,1); else, C = eye(P).*eps; end
        else
            if n>=1, C = diag(var(X,1,1)); else, C = eye(P).*eps; end
        end
        C = C + reg*eye(P);
    end
end

function plot_mahal_mean_sem(out, taxis)
% out: struct from mahal_from_trials
cols = lines(numel(out.pairs)); LW = 2.8; a = 0.22;
figure('Color','w'); hold on;
for k = 1:numel(out.pairs)
    shaded_band(taxis, out.D(:,k), 2*out.Dsem(:,k), cols(k,:), LW, a);
end
xlabel('Time'); ylabel('Mahalanobis distance');
legend({'Thr-LS', 'Thr-Ch', 'LS-Ch'},'Location','best'); box on; grid on;
title('Pairwise Mahalanobis distance (mean \pm 2 SEM)');
end

function shaded_band(t, y, ysem, color, lw, a)

t = t(:); y = y(:); ysem = ysem(:);
lo = y - ysem; hi = y + ysem;
patch([t; flipud(t)], [lo; flipud(hi)], color, ...
      'EdgeColor','none','FaceAlpha',a);
plot(t, y, 'Color', color, 'LineWidth', lw);
end

function make_marker_movies_meanfade(markerData2plot, varargin)
% Movies: mean marker positions with fading tail
%
% markerData2plot(bhv).data: [T x Ntrials x (2*M)]
% packed as [x1 y1 x2 y2 ...]
%
% Time options:
%   'Time'       : [T x 1] ms, 0 = onset
%   'MsPerFrame' : scalar ms per frame
%   'OnsetIdx'   : onset index (with MsPerFrame)
%
% Key options:
%   'TopFrac'      (default 0.20)
%   'MinTrials'    (default 10)
%   'NTrail'       (default 5)     % how many past frames are kept
%   'Colors'       (default [1 0 0; 0 0 1; 0 1 0])
%   'GestureNames' (default {'threat','lipsmack','chew'})
%   'InvertY'      (default true)
%   'FPS'          (default 24)
%   'FrameStep'    (default 1)
%   'OutDir'       (default pwd)

p = inputParser;
addParameter(p,'Time',[],@(x)isnumeric(x)&&isvector(x));
addParameter(p,'MsPerFrame',[],@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'OnsetIdx',[],@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'TopFrac',0.20);
addParameter(p,'MinTrials',10);
addParameter(p,'NTrail',5);
addParameter(p,'Colors',[1 0 0; 0 0 1; 0 1 0]);
addParameter(p,'GestureNames',{'threat','lipsmack','chew'});
addParameter(p,'InvertY',true);
addParameter(p,'FPS',24);
addParameter(p,'FrameStep',1);
addParameter(p,'OutDir',pwd);
parse(p,varargin{:});
opt = p.Results;

[T,~,P2] = size(markerData2plot(1).data);
M = P2/2;

% ----- time -----
if ~isempty(opt.Time)
    tms = opt.Time(:);
else
    tms = ((1:T)' - opt.OnsetIdx)*opt.MsPerFrame;
end
win = [-500 2000]; base = [-1000 -500]; post = [0 2000];
winIdx  = find(tms>=win(1)  & tms<=win(2));
baseIdx = find(tms>=base(1) & tms<=base(2));
postIdx = find(tms>=post(1) & tms<=post(2));

% ----- select top trials per gesture -----
idxTrials = cell(1,3);
for g=1:3
    Xg = markerData2plot(g).data; N=size(Xg,2);
    baseMean = squeeze(nanmean(Xg(baseIdx,:,:),1));
    postData = Xg(postIdx,:,:);
    magn=nan(N,1);
    for n=1:N
        Bn = baseMean(n,:);
        D  = squeeze(postData(:,n,:)) - Bn;
        d2 = sum(D.^2,2,'omitnan');
        magn(n)=sqrt(max(d2,[],'omitnan'));
    end
    k = max(opt.MinTrials, ceil(opt.TopFrac*N));
    [~,ord] = sort(magn,'descend','MissingPlacement','last');
    idxTrials{g} = sort(ord(1:k));
end

% ----- global axes -----
[xmin,xmax,ymin,ymax]=deal(+inf,-inf,+inf,-inf);
for g=1:3
    d=markerData2plot(g).data(winIdx, idxTrials{g},:);
    X=d(:,:,1:2:end); Y=d(:,:,2:2:end); if opt.InvertY, Y=-Y; end
    xmin=min(xmin,min(X(:))); xmax=max(xmax,max(X(:)));
    ymin=min(ymin,min(Y(:))); ymax=max(ymax,max(Y(:)));
end
dx=xmax-xmin; dy=ymax-ymin; pad=0.05;
axlim=[xmin-pad*dx xmax+pad*dx ymin-pad*dy ymax+pad*dy];

% ----- render movies -----
if ~exist(opt.OutDir,'dir'), mkdir(opt.OutDir); end

for g=1:3
    trials=idxTrials{g};
    d=markerData2plot(g).data(winIdx,trials,:);
    Tw=size(d,1);
    Xall=d(:,:,1:2:end); Yall=d(:,:,2:2:end);
    if opt.InvertY, Yall=-Yall; end
    Xmu=squeeze(nanmean(Xall,2)); Ymu=squeeze(nanmean(Yall,2));

    fig=figure('Color','w','Position',[100 100 800 800]);
    ax=axes('Parent',fig); hold(ax,'on'); axis(ax,'equal');
    xlim(ax,axlim(1:2)); ylim(ax,axlim(3:4));
    set(ax,'Visible','off');

    % text overlays
    xL=axlim(1); xR=axlim(2); yB=axlim(3); yT=axlim(4);
    xLabel=xL+0.02*(xR-xL); xClock=xL+0.68*(xR-xL); yText=yB+0.04*(yT-yB);
    text(ax,xLabel,yText,opt.GestureNames{g},'Color',opt.Colors(g,:), ...
        'FontSize',16,'FontWeight','bold','HorizontalAlignment','left');
    timeTxt=text(ax,xClock,yText,'','Color',[0 0 0], ...
        'FontSize',14,'FontWeight','bold','HorizontalAlignment','left');

    % scatter handles for trail points
    scatters=gobjects(opt.NTrail,M);
    for m=1:M
        for j=1:opt.NTrail
            scatters(j,m)=scatter(ax,nan,nan,80,opt.Colors(g,:),'filled', ...
                                  'MarkerFaceAlpha', (j/opt.NTrail)*0.8);
        end
    end

    vname=fullfile(opt.OutDir,sprintf('gesture_%d.mp4',g));
    vw=VideoWriter(vname,'MPEG-4'); vw.FrameRate=opt.FPS; open(vw);

    for tIdx=1:opt.FrameStep:Tw
        % update trail scatters
        for m=1:M
            for j=1:opt.NTrail
                t0=tIdx-(j-1);
                if t0>=1
                    set(scatters(j,m),'XData',Xmu(t0,m),'YData',Ymu(t0,m));
                else
                    set(scatters(j,m),'XData',nan,'YData',nan);
                end
            end
        end
        set(timeTxt,'String',sprintf('time = %d ms',round(tms(winIdx(tIdx)))));
        drawnow limitrate nocallbacks;
        writeVideo(vw,getframe(fig));
    end
    close(vw); close(fig);
    fprintf('Saved %s\n',vname);
end
end

function make_marker_movies(markerData2plot, varargin)

% MAKE_MARKER_MOVIES  Per-gesture movies of 2-D facial markers.
% - Selects top 20% biggest trials by displacement (baseline -500..0ms vs 0..2000ms)
% - Renders window -500..+2000 ms around onset
% - Hides axes/ticks
% - Adds bottom label (gesture name, same color) + running clock text per frame
%
% Required input:
%   markerData2plot(bhv).data : [T x Ntrials x (2*M)], packed [x1 y1 x2 y2 ...]
%
% Time options (choose one):
%   'Time'       : [T x 1] ms, 0 = movement onset
%   'MsPerFrame' : scalar (ms)
%   'OnsetIdx'   : scalar index of onset
%
% Other options:
%   'OutDir'        (default pwd)
%   'FPS'           (default 30)
%   'FrameStep'     (default 1)
%   'InvertY'       (default true)
%   'TopFrac'       (default 0.20)
%   'MinTrials'     (default 10)
%   'Colors'        (default [1 0 0; 0 0 1; 0 1 0])  % red, blue, green
%   'GestureNames'  (default {'threat','lipsmack','chew'})
%
% Outputs: gesture_1.mp4, gesture_2.mp4, gesture_3.mp4 in OutDir

p = inputParser;
addParameter(p,'Time',[],@(x)isnumeric(x)&&isvector(x));
addParameter(p,'MsPerFrame',[],@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'OnsetIdx',[],@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'OutDir',pwd,@(s)ischar(s)||isstring(s));
addParameter(p,'FPS',30,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'FrameStep',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'InvertY',true,@islogical);
addParameter(p,'TopFrac',0.20,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
addParameter(p,'MinTrials',10,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'Colors',[1 0 0; 0 0 1; 0 1 0],@(x)isnumeric(x)&&isequal(size(x),[3,3])); % R,B,G
addParameter(p,'GestureNames',{'threat','lipsmack','chew'},@(c)iscellstr(c)&&numel(c)==3);
parse(p,varargin{:});
opt = p.Results;

assert(numel(markerData2plot)==3,'Expect 1x3 structure array.');
[T, ~, P2] = size(markerData2plot(1).data);
assert(mod(P2,2)==0,'Positions dim must be even [x1 y1 ...].');
M = P2/2;

% --- Resolve time vector in ms (0 = onset) ---
if ~isempty(opt.Time)
    tms = opt.Time(:); assert(numel(tms)==T,'Time length must match T.');
else
    assert(~isempty(opt.MsPerFrame) && ~isempty(opt.OnsetIdx), ...
        'Provide Time OR MsPerFrame+OnsetIdx.');
    tms = ((1:T)' - opt.OnsetIdx) * opt.MsPerFrame;
end

% Windows (ms)
tWin  = [-500 2000];
tBase = [-1000 -500];
tPost = [0 2000];

winIdx  = find(tms>=tWin(1)  & tms<=tWin(2));
baseIdx = find(tms>=tBase(1) & tms<=tBase(2));
postIdx = find(tms>=tPost(1) & tms<=tPost(2));
assert(~isempty(winIdx)&&~isempty(baseIdx)&&~isempty(postIdx),'Window not in Time.');

% --- Select biggest trials per gesture (top 20% by max displacement) ---
idxTrials = cell(1,3);
for g=1:3
    Xg = markerData2plot(g).data;   % [T x N x 2M]
    N  = size(Xg,2);
    baseMean = squeeze(nanmean(Xg(baseIdx,:,:),1)); % [N x 2M]
    postData = Xg(postIdx,:,:);                      % [Tp x N x 2M]
    Tp = size(postData,1);
    magn = nan(N,1);
    for n=1:N
        Bn = baseMean(n,:);
        D  = squeeze(postData(:,n,:)) - Bn;   % [Tp x 2M]
        d2 = sum(D.^2, 2, 'omitnan');         % per-time squared norm across all markers
        magn(n) = sqrt(max(d2,[],'omitnan'));
    end
    k = max(opt.MinTrials, ceil(opt.TopFrac*N));
    [~,ord] = sort(magn,'descend','MissingPlacement','last');
    idxTrials{g} = sort(ord(1:k));
end

% --- Global axes limits across gestures (consistent framing) ---
[xmin,xmax,ymin,ymax] = deal(+inf,-inf,+inf,-inf);
for g=1:3
    d = markerData2plot(g).data(winIdx, idxTrials{g}, :);
    X = d(:,:,1:2:end); Y = d(:,:,2:2:end);
    if opt.InvertY, Y = -Y; end
    xmin = min(xmin, min(X(:),[],'omitnan')); xmax = max(xmax, max(X(:),[],'omitnan'));
    ymin = min(ymin, min(Y(:),[],'omitnan')); ymax = max(ymax, max(Y(:),[],'omitnan'));
end
pad = 0.05; dx = xmax-xmin; dy = ymax-ymin; if dx<=0, dx=1; end; if dy<=0, dy=1; end
axlim = [xmin - pad*dx, xmax + pad*dx, ymin - pad*dy, ymax + pad*dy];

% --- Output ---
if ~exist(opt.OutDir,'dir'), mkdir(opt.OutDir); end

% --- Render per gesture ---
for g=1:3
    trials = idxTrials{g};
    d  = markerData2plot(g).data(winIdx, trials, :); % [Tw x k x 2M]
    Tw = size(d,1); k = numel(trials);
    X  = d(:,:,1:2:end);   % [Tw x k x M]
    Y  = d(:,:,2:2:end);   % [Tw x k x M]
    if opt.InvertY, Y = -Y; end

    fig = figure('Color','w','Position',[100 100 800 800]);
    ax  = axes('Parent',fig); hold(ax,'on'); axis(ax,'equal');
    xlim(ax, axlim(1:2)); ylim(ax, axlim(3:4));
    set(ax,'Visible','off'); % hide axes/ticks

    % Bottom overlay positions (in data coords)
    xL = axlim(1); xR = axlim(2); yB = axlim(3); yT = axlim(4);
    xLabel = xL + 0.02*(xR-xL);   % left 2%
    xClock = xL + 0.68*(xR-xL);   % right-ish
    yText  = yB + 0.04*(yT-yB);   % 4% above bottom

    % Add static gesture label (same color)
    labelTxt = text(ax, xLabel, yText, opt.GestureNames{g}, ...
        'Color', opt.Colors(g,:), 'FontSize', 16, 'FontWeight', 'bold', ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
        'Interpreter','none');

    % Add time clock (updated every frame)
    timeTxt  = text(ax, xClock, yText, sprintf('time = %d ms', round(tms(winIdx(1)))), ...
        'Color', [0 0 0], 'FontSize', 14, 'FontWeight', 'bold', ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
        'Interpreter','none');

    % Precreate artists
    scatters = gobjects(M,1);
    xi0 = squeeze(X(1,:,:)); yi0 = squeeze(Y(1,:,:)); % [k x M]
    for m=1:M
        scatters(m) = scatter(ax, xi0(:,m), yi0(:,m), 36, opt.Colors(g,:), 'filled', ...
                              'MarkerFaceAlpha', 0.18, 'MarkerEdgeAlpha', 0.0);
    end

    vname = fullfile(opt.OutDir, sprintf('gesture_%d.mp4', g));
    vw = VideoWriter(vname, 'MPEG-4'); vw.FrameRate = opt.FPS; open(vw);

    for tIdx = 1:opt.FrameStep:Tw
        % Update points
        xi = squeeze(X(tIdx,:,:)); yi = squeeze(Y(tIdx,:,:)); % [k x M]
        for m=1:M
            set(scatters(m), 'XData', xi(:,m), 'YData', yi(:,m));
        end
        % Update time clock
        set(timeTxt, 'String', sprintf('time = %d ms', round(tms(winIdx(tIdx)))));
        drawnow limitrate nocallbacks;
        writeVideo(vw, getframe(fig));
    end

    close(vw); close(fig);
    fprintf('Saved %s (trials=%d, window=%d..%d ms)\n', vname, k, ...
        round(tms(winIdx(1))), round(tms(winIdx(end))));
end

fprintf('All movies done.\n');
end

function make_marker_movies_biggestTrial(markerData2plot, varargin)
% Movies: single-trial marker positions with fading tail
%
% markerData2plot(bhv).data: [T x Ntrials x (2*M)]
% packed as [x1 y1 x2 y2 ...]
%
% Time options:
%   'Time'       : [T x 1] ms, 0 = onset
%   'MsPerFrame' : scalar ms per frame
%   'OnsetIdx'   : onset index (with MsPerFrame)
%
% Plot options:
%   'NTrail'       (default 5)     % how many past frames are kept
%   'Colors'       (default [1 0 0; 0 0 1; 0 1 0])
%   'GestureNames' (default {'threat','lipsmack','chew'})
%   'InvertY'      (default true)
%   'FPS'          (default 24)
%   'FrameStep'    (default 1)
%   'OutDir'       (default pwd)

p = inputParser;
addParameter(p,'Time',[],@(x)isnumeric(x)&&isvector(x));
addParameter(p,'MsPerFrame',[],@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'OnsetIdx',[],@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'NTrail',5);
addParameter(p,'Colors',[1 0 0; 0 0 1; 0 1 0]);
addParameter(p,'GestureNames',{'threat','lipsmack','chew'});
addParameter(p,'InvertY',true);
addParameter(p,'FPS',24);
addParameter(p,'FrameStep',1);
addParameter(p,'OutDir',pwd);
parse(p,varargin{:});
opt = p.Results;

[T,~,P2] = size(markerData2plot(1).data);
M = P2/2;

% ----- time -----
if ~isempty(opt.Time)
    tms = opt.Time(:);
else
    tms = ((1:T)' - opt.OnsetIdx)*opt.MsPerFrame;
end
win  = [-500 2000];
base = [-1000 -500];
post = [0 2000];
winIdx  = find(tms>=win(1)  & tms<=win(2));
baseIdx = find(tms>=base(1) & tms<=base(2));
postIdx = find(tms>=post(1) & tms<=post(2));

% ----- select SINGLE fastest trial per gesture (speed-based) -----
% Requires MsPerFrame or Time to infer frame duration.
if isempty(opt.MsPerFrame)
    if ~isempty(opt.Time)
        % infer ms/frame from Time
        dt = mode(diff(opt.Time));
        MsPerFrame_eff = dt;
    else
        error('To select by speed, provide MsPerFrame or Time to infer frame duration.');
    end
else
    MsPerFrame_eff = opt.MsPerFrame;
end

idxTrial = nan(1,3);
for g = 1:3
    Xg = markerData2plot(g).data;           % [T x N x 2M]
    N  = size(Xg,2);
    if N==0
        idxTrial(g) = NaN; 
        continue;
    end

    % Use post-onset period for speed computation
    postData = Xg(postIdx,:,:);             % [Tpost x N x 2M]

    % speed per trial: mean Euclidean step size per frame / MsPerFrame
    spd = nan(N,1);
    for n = 1:N
        % positions over time for trial n in 2M-dim space
        P = squeeze(postData(:,n,:));       % [Tpost x 2M]
        if size(P,1) < 2 || all(~isfinite(P(:)))
            spd(n) = NaN; 
            continue;
        end
        dP = diff(P,1,1);                   % frame-to-frame steps [Tpost-1 x 2M]
        step = sqrt(sum(dP.^2,2));          % Euclidean step length per frame
        spd_per_frame = step ./ MsPerFrame_eff;  % units per ms (×1000 for per s)
        % Summary choice: mean speed (robust alternative: prctile(spd_per_frame,95))
        spd(n) = mean(spd_per_frame,'omitnan');
    end

    % pick the single fastest trial
    [~,ord] = sort(spd,'descend','MissingPlacement','last');
    idxTrial(g) = ord(2);
end


% ----- global axes limits from the chosen single trial in each gesture -----
[xmin,xmax,ymin,ymax] = deal(+inf,-inf,+inf,-inf);
for g = 1:3
    if isnan(idxTrial(g)), continue; end
    d  = markerData2plot(g).data(winIdx, idxTrial(g), :);   % [Twin x 1 x 2M]
    X  = d(:,:,1:2:end); Y = d(:,:,2:2:end);
    if opt.InvertY, Y = -Y; end
    xmin = min(xmin, min(X(:))); xmax = max(xmax, max(X(:)));
    ymin = min(ymin, min(Y(:))); ymax = max(ymax, max(Y(:)));
end
dx = xmax - xmin; dy = ymax - ymin; pad = 0.05;
if ~isfinite(dx) || ~isfinite(dy) || dx==0 || dy==0
    dx = 1; dy = 1; xmin = -0.5; xmax = 0.5; ymin = -0.5; ymax = 0.5;
end
axlim = [xmin - pad*dx, xmax + pad*dx, ymin - pad*dy, ymax + pad*dy];

% ----- render movies -----
if ~exist(opt.OutDir,'dir'), mkdir(opt.OutDir); end

for g = 1:3
    if isnan(idxTrial(g))
        warning('Gesture %d has no valid trials. Skipping.', g);
        continue;
    end

    d   = markerData2plot(g).data(winIdx, idxTrial(g), :);  % [Twin x 1 x 2M]
    Tw  = size(d,1);
    Xtr = squeeze(d(:,:,1:2:end));                          % [Twin x M]
    Ytr = squeeze(d(:,:,2:2:end));                          % [Twin x M]
    if opt.InvertY, Ytr = -Ytr; end

    fig = figure('Color','w','Position',[100 100 800 800]);
    ax  = axes('Parent',fig); hold(ax,'on'); axis(ax,'equal');
    xlim(ax,axlim(1:2)); ylim(ax,axlim(3:4));
    set(ax,'Visible','off');

    % text overlays
    xL=axlim(1); xR=axlim(2); yB=axlim(3); yT=axlim(4);
    xLabel = xL + 0.02*(xR-xL);
    xClock = xL + 0.68*(xR-xL);
    yText  = yB + 0.04*(yT-yB);
    text(ax,xLabel,yText,opt.GestureNames{g},'Color',opt.Colors(g,:), ...
        'FontSize',16,'FontWeight','bold','HorizontalAlignment','left');
    timeTxt = text(ax,xClock,yText,'','Color',[0 0 0], ...
        'FontSize',14,'FontWeight','bold','HorizontalAlignment','left');

    % trail scatter handles (single-trial) ------------------------------ %% CHANGED
    scatters = gobjects(opt.NTrail, M);
    for m = 1:M
        for j = 1:opt.NTrail
            scatters(j,m) = scatter(ax, nan, nan, 80, opt.Colors(g,:), ...
                'filled', 'MarkerFaceAlpha', (j/opt.NTrail)*0.8);
        end
    end

    vname = fullfile(opt.OutDir, sprintf('gesture_%d_biggestTrial.mp4', g));
    vw = VideoWriter(vname,'MPEG-4'); vw.FrameRate = opt.FPS; open(vw);

    for tIdx = 1:opt.FrameStep:Tw
        % update trail points from the SINGLE trial trajectory  ---------- %% CHANGED
        for m = 1:M
            for j = 1:opt.NTrail
                t0 = tIdx - (j-1);
                if t0 >= 1
                    set(scatters(j,m), 'XData', Xtr(t0,m), 'YData', Ytr(t0,m));
                else
                    set(scatters(j,m), 'XData', nan, 'YData', nan);
                end
            end
        end
        set(timeTxt,'String', sprintf('time = %d ms', round(tms(winIdx(tIdx)))));
        drawnow limitrate nocallbacks;
        writeVideo(vw, getframe(fig));
    end

    close(vw); close(fig);
    fprintf('Saved %s\n', vname);
end
end
