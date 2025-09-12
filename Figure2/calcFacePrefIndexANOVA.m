
%% run 2 way anova for effect of cortical region, gesture type, on mean
%% FR & run post-hoc tests
%% these were non-significant (reported but not plotted) in main text 

%% run 1 way anova for effect of cortical region on FEP 
%% and run post-hoc tests -- plotted and reported 

clear all; 
set(0,'defaultAxesFontSize', 16); % bc im blind 
subjects = {'barney','thor'};

pathToMatFiles = '/Users/geena/Dropbox/PhD/SUAinfo/';

%% get barney's data 
barneyRawData1 = load(fullfile(pathToMatFiles, 'Barney_210704/Data4Analysis/210704_Barney_faceExpPref_V3.mat'));
barneyRawData2 = load(fullfile(pathToMatFiles, 'Barney_210706/Data4Analysis/210706_Barney_faceExpPref_V3.mat'));
barneyRawData3 = load(fullfile(pathToMatFiles, 'Barney_210805/Data4Analysis/210805_Barney_faceExpPref_V3.mat'));

barneyRawData.FEP = [barneyRawData1.faceExpPref, barneyRawData2.faceExpPref, barneyRawData3.faceExpPref]';
barneyRawData.meanFRs = cat(1, barneyRawData1.Ris, barneyRawData2.Ris, barneyRawData3.Ris);
barneyRawData.spikeLabels2plot = [barneyRawData1.spikeLabels2plot; barneyRawData2.spikeLabels2plot;  barneyRawData3.spikeLabels2plot]';
barneyRawData.cortRegion = [barneyRawData1.cortRegion barneyRawData2.cortRegion barneyRawData3.cortRegion]';

%% get thor's data 
thorRawData1 = load(fullfile(pathToMatFiles, 'Thor_171010/Data4Analysis/171010_Thor_faceExpPref_V3.mat'));
thorRawData2 = load(fullfile(pathToMatFiles, 'Thor_171027/Data4Analysis/171027_Thor_faceExpPref_V3.mat'));
thorRawData3 = load(fullfile(pathToMatFiles, 'Thor_171005/Data4Analysis/171005_Thor_faceExpPref_V3.mat'));
thorRawData4 = load(fullfile(pathToMatFiles, 'Thor_171128/Data4Analysis/171128_Thor_faceExpPref_V3.mat'));

thorRawData.FEP = [thorRawData1.faceExpPref, thorRawData2.faceExpPref, thorRawData3.faceExpPref, thorRawData4.faceExpPref]';
thorRawData.meanFRs = cat(1, thorRawData1.Ris, thorRawData2.Ris, thorRawData3.Ris, thorRawData4.Ris);
thorRawData.spikeLabels2plot = [thorRawData1.spikeLabels2plot; thorRawData2.spikeLabels2plot;  thorRawData3.spikeLabels2plot; thorRawData4.spikeLabels2plot]';
thorRawData.cortRegion = [thorRawData1.cortRegion thorRawData2.cortRegion thorRawData3.cortRegion thorRawData4.cortRegion]';

%% get the cortical regions per cell 
cortRegionOut = [];
FEPout = [];
meanFRsOut = [];

for ss = 1:length(subjects)
    if strcmpi(subjects{ss},'Barney')
        data = barneyRawData;
    elseif strcmpi(subjects{ss},'Thor')
        data = thorRawData;
    end
    
    % pool faceExpPrefs indices 
    FEP = data.FEP;
    cortRegion = data.cortRegion; 
    meanFRs = data.meanFRs;

    cortRegionOut = [cortRegionOut; cortRegion];
    FEPout = [FEPout ; FEP];
    meanFRsOut = [meanFRsOut; meanFRs];

    clear cortRegion; clear FEP; clear meanFRs;
end

% run 1 way anova for effect of cortical region on FEP 
% and run post-hoc tests 
[P_oneway,ANOVATAB_oneway,STATS_oneway] = anova1(FEPout, cortRegionOut,'off');
[c_oneway,m_oneway] = multcompare(STATS_oneway);

% run 2 way anova for effect of cortical region, expresison type, on mean
% FR & run post-hoc tests

inputFR = [meanFRsOut(:,1) ; meanFRsOut(:,2); meanFRsOut(:,3)];
T    = cell(1, size(meanFRsOut,1)); T(:) = {'Threat'}; L = cell(1, size(meanFRsOut,1)); L(:) = {'Lipsmack'}; C = cell(1, size(meanFRsOut,1)); C(:) = {'Chew'};
inputExp = [T'; L'; C'];
inputCortRegion = [cortRegionOut ; cortRegionOut; cortRegionOut];
[P_twoway,TABLE_twoway,STATS_twoway] = anovan(inputFR,{inputExp inputCortRegion},'model','interaction','varnames',{'inputExp','inputCortRegion'});
[c_twoway,m_twoway, ~, names] = multcompare(STATS_twoway,'Dimension', [1 2]);

% plot it 
factors = STATS_twoway.varnames;
groupNames = STATS_twoway.grpnames;
cmap = [0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	];
start = 1; stop = 3;
for rr = 1:length(groupNames{2,1}) % per region
    
    means = m_twoway(start:stop,1);
    errors = m_twoway(start:stop,2);
    
    b(rr) = bar(start:1:stop, means); 
    b(rr).FaceColor = cmap(rr,:);
    hold on

    er = errorbar(start:stop, means,errors,'HandleVisibility','off');
    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 4;

    start = start +3; stop = stop + 3;
end
hold off

tickLabels = cell(1, length(m_twoway)); tickLabels(:) = {'Threat','Lipsmack','Chew','Threat','Lipsmack','Chew','Threat','Lipsmack','Chew','Threat','Lipsmack','Chew'};
set(gca,'XTick',1:length(m_twoway), 'XTickLabel',tickLabels);
legend({'S1','M1','PMv','M3'});

title('MeanFR by Region & GestureType')
ylabel('mean FR');

%% plot again just colored differently
cmap2 = [1 0 0
    0 0 1
    0 1 0];
start = 1; skip = 3;
meansAll = [];
errorsAll = [];
for exp = 1:length(groupNames{1,1}) % per region
    max = length(m_twoway);

    means = m_twoway(start:skip:max,1);
    errors = m_twoway(start:skip:max,2);

    meansAll = [meansAll means];
    errorsAll = [errorsAll errors];
    start = start +1;

end
%%
figure; 
b = bar(meansAll); hold on;

ngroups = size(meansAll, 1);
nbars = size(meansAll, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, meansAll(:,i), errorsAll(:,i),'HandleVisibility','off');
    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 4;
end
hold off

b(1).FaceColor = cmap2(1,:);b(2).FaceColor = cmap2(2,:);b(3).FaceColor = cmap2(3,:);

tickLabels = cell(1, 4); tickLabels(:) = {'S1','M1','PMv','M3'};
set(gca,'XTick',1:length(m_twoway), 'XTickLabel',tickLabels);
legend({'Threat','Lipsmack','Chew'});

title('MeanFR by Region & GestureType')
ylabel('mean FR');

%% === APA-style reporting for two-way ANOVA + Tukey post-hoc (no local functions) ===
%% === APA-style reporting for two-way ANOVA + Tukey post-hoc ===

% Identify columns in the ANOVA table
hdr   = string(TABLE_twoway(1,:));
colSS = find(contains(hdr, ["Sum","SS"], 'IgnoreCase', true), 1);
colDF = find(contains(hdr, ["d.f","df"], 'IgnoreCase', true), 1);
colF  = find(strcmpi(hdr, "F"), 1);
colP  = find(contains(hdr, ["Prob>F","p"], 'IgnoreCase', true), 1);

% Row indices
rowNames = string(TABLE_twoway(2:end,1));                   % names without header
getRow   = @(lab) find(strcmp(rowNames, string(lab)), 1) + 1;  % +1 for header row

% Factor labels (match your anovan varnames)
facA   = 'inputExp';
facB   = 'inputCortRegion';
facInt = [facA '*' facB];

rA   = getRow(facA);
rB   = getRow(facB);
rInt = getRow(facInt);
rErr = find(contains(rowNames, ["Error","Residual"], 'IgnoreCase', true), 1) + 1;

% Friendly print labels
labelA   = 'Expression';
labelB   = 'Cortical region';
labelInt = 'Expression × Cortical region';

disp('Two-way ANOVA (APA-style):');

if ~isempty(rA) && ~isempty(rErr)
    df1 = TABLE_twoway{rA, colDF}; df2 = TABLE_twoway{rErr, colDF};
    Fv  = TABLE_twoway{rA, colF};  p   = TABLE_twoway{rA, colP};
    fprintf('%s: F(%d, %d) = %.2f, %s\n', labelA, df1, df2, Fv, formatP(p));
end

if ~isempty(rB) && ~isempty(rErr)
    df1 = TABLE_twoway{rB, colDF}; df2 = TABLE_twoway{rErr, colDF};
    Fv  = TABLE_twoway{rB, colF};  p   = TABLE_twoway{rB, colP};
    fprintf('%s: F(%d, %d) = %.2f, %s\n', labelB, df1, df2, Fv, formatP(p));
end

if ~isempty(rInt) && ~isempty(rErr)
    df1 = TABLE_twoway{rInt, colDF}; df2 = TABLE_twoway{rErr, colDF};
    Fv  = TABLE_twoway{rInt, colF};  p   = TABLE_twoway{rInt, colP};
    fprintf('%s: F(%d, %d) = %.2f, %s\n', labelInt, df1, df2, Fv, formatP(p));
end

% Optional: partial eta-squared for each effect
if ~isempty(colSS) && ~isempty(rErr)
    SSerr = TABLE_twoway{rErr, colSS};
    if ~isempty(rA)
        eta = TABLE_twoway{rA,colSS} / (TABLE_twoway{rA,colSS} + SSerr);
        fprintf('  eta_p^2 (%s) = %.3f\n', labelA, eta);
    end
    if ~isempty(rB)
        eta = TABLE_twoway{rB,colSS} / (TABLE_twoway{rB,colSS} + SSerr);
        fprintf('  eta_p^2 (%s) = %.3f\n', labelB, eta);
    end
    if ~isempty(rInt)
        eta = TABLE_twoway{rInt,colSS} / (TABLE_twoway{rInt,colSS} + SSerr);
        fprintf('  eta_p^2 (%s) = %.3f\n', labelInt, eta);
    end
end

%% Post-hoc (Tukey–Kramer) comparisons from multcompare
% c_twoway rows: [i j lower meanDiff upper pAdj]
% names: cell array of group labels for each (Expression × Region) cell
disp('Post-hoc (Tukey–Kramer) comparisons:');
for r = 1:size(c_twoway,1)
    i  = c_twoway(r,1);
    j  = c_twoway(r,2);
    lo = c_twoway(r,3);
    md = c_twoway(r,4);
    hi = c_twoway(r,5);
    p  = c_twoway(r,6);

    % Ensure names{i}/names{j} printable as char
    ni = char(string(names{i}));
    nj = char(string(names{j}));

    fprintf('%s vs %s: ΔM = %.3f, 95%% CI [%.3f, %.3f], %s (Tukey–Kramer)\n', ...
        ni, nj, md, lo, hi, formatP(p));
end

%% (Optional) Write the same lines to a UTF-8 text file
fid = fopen('anova_tukey_report.txt','w','n','UTF-8');
if fid > 0
    if ~isempty(rA) && ~isempty(rErr)
        df1 = TABLE_twoway{rA, colDF}; df2 = TABLE_twoway{rErr, colDF};
        Fv  = TABLE_twoway{rA, colF};  p   = TABLE_twoway{rA, colP};
        fprintf(fid,'%s: F(%d, %d) = %.2f, %s\n', labelA, df1, df2, Fv, formatP(p));
    end
    if ~isempty(rB) && ~isempty(rErr)
        df1 = TABLE_twoway{rB, colDF}; df2 = TABLE_twoway{rErr, colDF};
        Fv  = TABLE_twoway{rB, colF};  p   = TABLE_twoway{rB, colP};
        fprintf(fid,'%s: F(%d, %d) = %.2f, %s\n', labelB, df1, df2, Fv, formatP(p));
    end
    if ~isempty(rInt) && ~isempty(rErr)
        df1 = TABLE_twoway{rInt, colDF}; df2 = TABLE_twoway{rErr, colDF};
        Fv  = TABLE_twoway{rInt, colF};  p   = TABLE_twoway{rInt, colP};
        fprintf(fid,'%s: F(%d, %d) = %.2f, %s\n', labelInt, df1, df2, Fv, formatP(p));
    end
    fprintf(fid,'Post-hoc (Tukey–Kramer) comparisons:\n');
    for r = 1:size(c_twoway,1)
        i  = c_twoway(r,1);
        j  = c_twoway(r,2);
        lo = c_twoway(r,3);
        md = c_twoway(r,4);
        hi = c_twoway(r,5);
        p  = c_twoway(r,6);
        ni = char(string(names{i}));
        nj = char(string(names{j}));
        fprintf(fid,'%s vs %s: ΔM = %.3f, 95%%%% CI [%.3f, %.3f], %s (Tukey–Kramer)\n', ...
            ni, nj, md, lo, hi, formatP(p));
    end
    fclose(fid);
end

%% ===== Helper placed at the end of the script =====
function pText = formatP(p)
%FORMATP  APA-style p-value formatting (no leading zero; < .001 if small)
    if ~isscalar(p) || ~isfinite(p)
        pText = 'p = NaN';
        return
    end
    if p < 0.001
        pText = 'p < .001';
    else
        pText = sprintf('p = %.3f', p);
        % remove leading zero
        pText = regexprep(pText, 'p = 0\.', 'p = .');
    end
end
