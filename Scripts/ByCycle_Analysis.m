%% ByCycle analysis: cycle-by-cycle features of OMEGA oscillatory episodes

%% Paths
p.bycycle = 'Z:\Enrique\Descriptives\Results\ByCycle'; % cycle features from ByCycle_OMEGA
p.episodes = 'Z:\Enrique\Descriptives\Results\Episodes';
p.figures = 'Z:\Enrique\Descriptives\Figures\ByCycle';
p.suppdata = 'Z:\Enrique\Descriptives\Supplementary Data'; % source data
p.files = 'Z:\Enrique\Descriptives\Files';

% Toolbox paths
addpath('Z:\Toolbox\fieldtrip-20230118') % Path to Fieldtrip
ft_defaults
addpath('Z:\Toolbox\slanCM\slanCM') % colors to plot
addpath(p.files)

%% OMEGA
% Participant list
sub = [1 2 3 4 5 6 7 8 9 11 12 14 15 16 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 37 39 40 41 42 44 45 46 47 48 49 50 51 52 55 56 57 58 59 60 61 62 63 64 65 67 68 69 70 71 72 73 74 75 76 77 78 79 80 84 85 87 88 89 90 91 92 94 95 96 97 98 99 101 102 103 104 105 106 134 145 146 148 149 150 151 152 154 155 156 157 158 159 160 161 165 166 167 168 169 170 171 175 176 177 179 181 184 185 195 197 200 207 208 210 212]';

Nsub  = length(sub);

load source_template_10mm_1925.mat
Nvoxin = length(find(source.inverse.inside));

% Frequency bands
freqnames = {'delta', 'theta', 'alpha', 'lowbeta', 'highbeta'};
Nbands = length(freqnames);

% Cycle features 
load([p.bycycle '\delta\cycle_features_s1.mat'])
features = fieldnames(cycle_features{1});
paramlist = {'band_amp', 'time_ptsym', 'time_rdsym', 'period', 'volt_amp'};
[~, paramidx] = ismember(paramlist, features);
[~, idx_rdsym] = ismember('time_rdsym', features);
[~, idx_ptsym] = ismember('time_ptsym', features);

% ROIs
load(fullfile(p.files, 'aal_voxel_label_10mm.mat'))
Nroi = length(aal_label_reduc);

% Permutations
nperm_t    = 1000;
nperm_corr = 5000;

% Figure settings
piecol  = jet_omega_mod;
piecols = [piecol(1,:); piecol(20,:); piecol(40,:); piecol(55,:); piecol(64,:)];
Ncols = 8;
Nrows = 5;

cmap = slanCM('vik', 256); 
cmap(100:156, :) = [];
cmap_no_white = interp1(linspace(0, 1, size(cmap,1)), cmap, linspace(0, 1, 256));
cmapC  = slanCM('vik'); % Figure 7
nC     = size(cmapC, 1);
grayNS = [0.8 0.8 0.8];

%% Avgs: median of each cycle feature per participant and voxel
median_features = NaN(Nsub, Nvoxin, length(features));
rng(1);
for fb = 1:Nbands
    frqband = freqnames{fb};
    display(frqband)
    for s = 1:Nsub
        load([p.bycycle '\' frqband '\cycle_features_s' num2str(s) '.mat'])
        for vx = 1:length(cycle_features)
            median_features(s,vx,:) = structfun(@(x) median(double(x)), cycle_features{vx});

            median_features(s,vx,idx_rdsym) = NaN;
            median_features(s,vx,idx_ptsym) = NaN;
            lt = double(cycle_features{vx}.sample_last_trough(:));
            if isempty(lt), continue; end
            nt    = double(cycle_features{vx}.sample_next_trough(:));
            epId  = cumsum([true; lt(2:end) ~= nt(1:end-1)]); 
            valid = find(accumarray(epId, 1) >= 3);
            if numel(valid) < 3, continue; end
            pt = double(cycle_features{vx}.time_ptsym(:));
            rd = double(cycle_features{vx}.time_rdsym(:));
            v  = nan(numel(valid), 2);
            for k = 1:numel(valid)
                idx = find(epId == valid(k));
                st  = randi(numel(idx) - 2);
                w   = idx(st:st + 2);
                v(k, :) = [abs(median(rd(w)) - 0.5) * 2, abs(median(pt(w)) - 0.5) * 2];
            end
            median_features(s,vx,idx_rdsym) = median(v(:, 1));
            median_features(s,vx,idx_ptsym) = median(v(:, 2));
        end
    end
    save(fullfile(p.bycycle, ['median_features_' freqnames{fb} '.mat']), 'median_features');
end

%% Statistics: ttest
fig6_paramlist = paramlist;
fig6_t    = nan(Nvoxin, Nbands, length(paramlist));
fig6_tsig = nan(Nvoxin, Nbands, length(paramlist));
fig6_p    = nan(Nvoxin, Nbands, length(paramlist));

for fb = 1:Nbands
    load([p.bycycle '\median_features_' freqnames{fb}])
    for param = 1:length(paramlist)
        % Load avg data
        tmpdata = squeeze(median_features(:,:,paramidx(param)));
        zData = normalize(tmpdata, 2, 'zscore');

        % Ttest
        [~, ~, ~, stats] = arrayfun(@(v) ttest(zData(:,v), 0), 1:size(tmpdata, 2),'UniformOutput',false);
        tVals = cellfun(@(s) s.tstat, stats);

        % Permutation test
        rng(1); 
        maxTval = zeros(nperm_t, 1);
        for pp = 1:nperm_t
            permsign = (rand(Nsub, 1) > 0.5)*2 - 1;
            permdata = zData .* permsign;
            [~,~,~, permstats] = arrayfun(@(v) ttest(permdata(:,v), 0), 1:size(permdata, 2),'UniformOutput',false);
            permtVals = cellfun(@(s) s.tstat, permstats);
            maxTval(pp) = max(abs(permtVals));
        end

        permTcorr = arrayfun(@(t) (sum(maxTval >= t) + 1) / (nperm_t + 1), abs(tVals));
        tValspermcor = tVals;
        tValspermcor(find(permTcorr > 0.05)) = 0;

        fig6_t(:, fb, param)    = tVals(:);
        fig6_tsig(:, fb, param) = tValspermcor(:);
        fig6_p(:, fb, param)    = permTcorr(:);

        % Figure 6: Cycle-by-cycle characterization
        cfg = [];
        cfg.colmap = cmap_no_white;
        cfg.colim = [-10 10];
        cfg.interp = 'linear';
        sBOSC_sourcefig(tValspermcor, cfg)
        exportgraphics(gcf, fullfile(p.figures, ['Figure6_' paramlist{param} '_' freqnames{fb} '.png']), 'Resolution', 300);
    end
end

%% ROIs: symmetry features
rdsym_roi = nan(Nsub, Nroi, Nbands);
ptsym_roi = nan(Nsub, Nroi, Nbands);
for fb = 1:Nbands
    load([p.bycycle '\median_features_' freqnames{fb}])
    rdsym_mat = squeeze(median_features(:, :, idx_rdsym));
    ptsym_mat = squeeze(median_features(:, :, idx_ptsym));

    for roi = 1:Nroi
        inds = find(label_inside_aal_reduc == roi);
        voxs = voxel_inside_aal(inds);
        rdsym_roi(:, roi, fb) = mean(rdsym_mat(:, voxs), 2, 'omitnan');
        ptsym_roi(:, roi, fb) = mean(ptsym_mat(:, voxs), 2, 'omitnan');
    end
end

%% Supplementary Figure 2 and 3
aal_short = {'PreCG', 'SFG', 'SFG (orb)', 'MFG', 'MFG (orb)', 'IFG (oper)', 'IFG (tri)', 'IFG (orb)', 'Rolandic', 'SMA', 'Olfactory', 'SFG (med)', 'SFG (med orb)', 'Rectus', 'Insula', 'ACG', 'MCG', 'PCG', 'Hippoc', 'ParaHippoc', 'Calcarine', 'Cuneus', 'Lingual', 'SOG', 'MOG', 'IOG', 'Fusiform', 'PostCG', 'SPG', 'IPG', 'Supramarg', 'Angular', 'Precuneus', 'ParaCL', 'Heschl', 'STG', 'Temp pole (sup)', 'MTG', 'Temp pole (mid)', 'ITG'};
aal_label_sort = {'PreCG', 'SFG', 'SFG (orb)', 'MFG', 'MFG (orb)', 'IFG (oper)', 'IFG (orb)',  'IFG (tri)', 'Rolandic', 'PostCG', 'SPG', 'IPG', 'Supramarg', 'Angular', 'PCG', 'ParaCL', 'Hippoc', 'ParaHippoc', 'Temp pole (mid)', 'Temp pole (sup)', 'STG', 'MTG', 'Heschl', 'ITG',  'Precuneus', 'Calcarine', 'Cuneus', 'Lingual', 'SOG', 'MOG', 'IOG', 'Fusiform', 'SMA', 'Olfactory', 'SFG (med)', 'SFG (med orb)', 'Rectus', 'Insula', 'ACG', 'MCG'};
[~, sorted_aal] = ismember(aal_label_sort, aal_short);
aal_short = aal_short(sorted_aal);
rdsym_roi = rdsym_roi(:,sorted_aal,:);
ptsym_roi = ptsym_roi(:,sorted_aal,:);

sym_roi   = {rdsym_roi, ptsym_roi};
sym_names = {'rdsym', 'ptsym'};
sym_xlim  = {[0.10 0.25], [0.08 0.20]};
for fb = 1:Nbands
    for sy = 1:length(sym_names)
        figure('Position', [100 100 Ncols*180 Nrows*180], 'Color', 'w');
        tiledlayout(Nrows, Ncols, 'TileSpacing', 'compact', 'Padding', 'tight');

        for roi = 1:Nroi
            x = sym_roi{sy}(:, roi, fb);
            valid = ~isnan(x);

            nexttile;
            histogram(x(valid), 10, 'FaceColor', piecols(fb,:), ...
                'EdgeColor', 'none', 'FaceAlpha', 0.7, 'Normalization', 'probability');
            hold on;
            xline(median(x(valid), 'omitnan'), 'k-', 'LineWidth', 1.5);
            xline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
            title(aal_short{roi}, 'Interpreter', 'none', 'FontSize', 10);
            set(gca, 'FontSize', 2.5, 'Box', 'off', 'TickDir', 'out');
            xlim(sym_xlim{sy})
            hold off;
        end

        fig = gcf;
        fig.Units = 'centimeters';
        fig.Position(3) = 9;
        fig.Position(4) = 5;

        exportgraphics(fig, fullfile(p.figures, [sym_names{sy} '_dist_roi_' freqnames{fb} '.png']), 'Resolution', 300);
        close;
    end
end

%% Correlations: within-subject spatial correlation between features 
featLabels = {'Duration', 'Power', 'Amplitude', 'Period', 'rdsym', 'ptsym'};
nFeat = numel(featLabels);

load(fullfile(p.episodes, 'dur_subj_vox.mat'))
load(fullfile(p.episodes, 'pow_subj_vox.mat'))

nSubj = size(dur_subj_vox, 1);

R_sub   = nan(nSubj, nFeat, nFeat, Nbands);   % correlation matrix
R_group = nan(nFeat, nFeat, Nbands);

for fb = 1:Nbands
    load([p.bycycle '\median_features_' freqnames{fb}])
    featuresdata = cat(3, ...
        dur_subj_vox(:, :, fb), ...                 % duration
        pow_subj_vox(:, :, fb), ...                 % power
        median_features(:, :, paramidx(1)), ...     % amplitude
        median_features(:, :, paramidx(4)), ...     % period
        median_features(:, :, paramidx(3)), ...     % rdsym
        median_features(:, :, paramidx(2)));        % ptsym

    % Correlate every pair of maps across voxels
    for s = 1:nSubj
        X = squeeze(featuresdata(s, :, :));
        R_sub(s, :, :, fb) = corr(X, 'type', 'Spearman', 'rows', 'pairwise');
    end

    % Fisher z and group average
    Z = atanh(R_sub(:, :, :, fb));
    Z(isinf(Z)) = NaN;
    Z2 = reshape(Z, nSubj, nFeat * nFeat);
    R_group(:, :, fb) = reshape(tanh(mean(Z2, 1, 'omitnan')), nFeat, nFeat);
    R_group(:, :, fb) = R_group(:, :, fb) - diag(diag(R_group(:, :, fb))) + eye(nFeat);
end

%  t-tests permutation
rng(1);
upmask = triu(true(nFeat), 1);
nU     = sum(upmask(:));
Zall   = [];
for fb = 1:Nbands
    Zfb  = atanh(R_sub(:, :, :, fb));
    Zfb(isinf(Zfb)) = NaN;
    Zfb  = reshape(Zfb, nSubj, nFeat * nFeat);
    Zall = [Zall, Zfb(:, upmask(:))];
end
nObs = sum(~isnan(Zall), 1);
tObs = mean(Zall, 1, 'omitnan') ./ (std(Zall, 0, 1, 'omitnan') ./ sqrt(nObs));
maxT = zeros(nperm_corr, 1);
for pp = 1:nperm_corr
    sgn      = (rand(nSubj, 1) > 0.5) * 2 - 1;
    Zp       = Zall .* sgn;
    tP       = mean(Zp, 1, 'omitnan') ./ (std(Zp, 0, 1, 'omitnan') ./ sqrt(nObs));
    maxT(pp) = max(abs(tP));
end
pPermAll = ((sum(maxT >= abs(tObs), 1) + 1) / (nperm_corr + 1))';

P_perm = nan(nFeat, nFeat, Nbands);
k = 0;
for fb = 1:Nbands
    tmpP = nan(nFeat); tmpP(upmask) = pPermAll(k + 1 : k + nU);
    P_perm(:, :, fb) = triu(tmpP, 1) + triu(tmpP, 1)';
    k = k + nU;
end

%% Figure 7: correlations between features (by band)
for fb = 1:Nbands
    fig = figure('Color', 'w', 'Position', [200 200 760 680], 'Renderer', 'painters');
    imagesc(R_group(:, :, fb), [-1 1]); colormap(cmapC); axis square;
    set(gca, 'XTick', 1:nFeat, 'XTickLabel', featLabels, 'YTick', 1:nFeat, ...
             'YTickLabel', featLabels, 'XTickLabelRotation', 90, 'FontSize', 13)
    axis off
    for i = 1:nFeat
        for j = 1:nFeat
            if i == j, continue; end
            r = R_group(i, j, fb);
            if P_perm(i, j, fb) < 0.05
                bg  = cmapC(round((min(max(r, -1), 1) + 1) / 2 * (nC - 1)) + 1, :);
                lum = 0.299 * bg(1) + 0.587 * bg(2) + 0.114 * bg(3);
                if lum < 0.5, fg = [1 1 1]; else, fg = [0 0 0]; end
                text(j, i, sprintf('%.2f', r), 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', fg)
            else
                rectangle('Position', [j - 0.5, i - 0.5, 1, 1], 'FaceColor', grayNS, 'EdgeColor', 'none')
            end
        end
    end
    exportgraphics(fig, fullfile(p.figures, ['feature_spatial_correlations_' freqnames{fb} '.png']), 'Resolution', 300)
    print(fig, fullfile(p.figures, ['feature_spatial_correlations_' freqnames{fb} '.svg']), '-dsvg')
end

%% Single correlation matrix across the 5 bands
Zb = atanh(R_sub);
Zb(isinf(Zb)) = NaN;
Zsub = mean(Zb, 4, 'omitnan');   % average
Z2 = reshape(Zsub, nSubj, nFeat * nFeat);
R_all = reshape(tanh(mean(Z2, 1, 'omitnan')), nFeat, nFeat);
R_all = R_all - diag(diag(R_all)) + eye(nFeat);

% t-tests  permutation
rng(1);
Zall  = Z2(:, upmask(:));
nObs  = sum(~isnan(Zall), 1);
tObs  = mean(Zall, 1, 'omitnan') ./ (std(Zall, 0, 1, 'omitnan') ./ sqrt(nObs));
maxT  = zeros(nperm_corr, 1);
for pp = 1:nperm_corr
    sgn      = (rand(nSubj, 1) > 0.5) * 2 - 1;
    Zp       = Zall .* sgn;
    tP       = mean(Zp, 1, 'omitnan') ./ (std(Zp, 0, 1, 'omitnan') ./ sqrt(nObs));
    maxT(pp) = max(abs(tP));
end
tmp = nan(nFeat); tmp(upmask) = (sum(maxT >= abs(tObs), 1) + 1) / (nperm_corr + 1);
P_all_perm = triu(tmp, 1) + triu(tmp, 1)';


%% Figure 7: correlations between features (all bands)
fig = figure('Color', 'w', 'Position', [200 200 760 680], 'Renderer', 'painters');
imagesc(R_all, [-1 1]); colormap(cmapC); axis square;
set(gca, 'XTick', 1:nFeat, 'XTickLabel', featLabels, 'YTick', 1:nFeat, ...
         'YTickLabel', featLabels, 'XTickLabelRotation', 90, 'FontSize', 13)
axis off
for i = 1:nFeat
    for j = 1:nFeat
        if i == j, continue; end
        r = R_all(i, j);
        if P_all_perm(i, j) < 0.05
            bg  = cmapC(round((min(max(r, -1), 1) + 1) / 2 * (nC - 1)) + 1, :);
            lum = 0.299 * bg(1) + 0.587 * bg(2) + 0.114 * bg(3);
            if lum < 0.5, fg = [1 1 1]; else, fg = [0 0 0]; end
            text(j, i, sprintf('%.2f', r), 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', fg)
        else
            rectangle('Position', [j - 0.5, i - 0.5, 1, 1], 'FaceColor', grayNS, 'EdgeColor', 'none')
        end
    end
end
exportgraphics(fig, fullfile(p.figures, 'feature_spatial_correlations_allbands.png'), 'Resolution', 300)
print(fig, fullfile(p.figures, 'feature_spatial_correlations_allbands.svg'), '-dsvg')

%% Save Source Data for Figures
save(fullfile(p.suppdata, 'SourceData_Fig6_7.mat'), 'fig6_paramlist', 'fig6_t', 'fig6_tsig', 'fig6_p', ...
     'featLabels', 'R_group', 'P_perm', 'R_all', 'P_all_perm');
