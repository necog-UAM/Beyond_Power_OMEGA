%% This script 

p.data = 'Z:\OMEGA\OMEGA_data';
p.results = 'Z:\Enrique\Waveform\Results';
p.figures = 'Z:\Enrique\Waveform\Figures';

% Toolbox paths
addpath('Z:\Toolbox\fieldtrip-20230118') % Path to Fieldtrip
addpath('Z:\Toolbox\NECOG')
ft_defaults

addpath(genpath('Z:\Toolbox\sourceBOSC'))

%% OMEGA 
% Participant and session list
sub = [1 2 3 4 5 6 7 8 9 11 12 14 15 16 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 37 39 40 41 42 44 45 46 47 48 49 50 51 52 55 56 57 58 59 60 61 62 63 64 65 67 68 69 70 71 72 73 74 75 76 77 78 79 80 84 85 87 88 89 90 91 92 94 95 96 97 98 99 101 102 103 104 105 106 134 145 146 148 149 150 151 152 154 155 156 157 158 159 160 161 165 166 167 168 169 170 171 175 176 177 179 181 184 185 195 197 200 207 208 210 212]';
ses = [1 1 1 1 1 1 1 2 1  2  1  2  1  1  1  2  3  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  3  1  2  1  2  1  2  1  1  2  1  1  1  1  1  1  1  1  2  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 ]';

subs = {} ; sess = {};
for s = 1:length(sub)
    if sub(s) < 10
        subs{s} = ['000' num2str(sub(s))];
    elseif sub(s) >= 10 & sub(s) < 100
        subs{s} = ['00' num2str(sub(s))];
    else
        subs{s} = ['0' num2str(sub(s))];
    end
    sess{s} = ['000' num2str(ses(s))];
end

Nsub  = length(sub);

load source_template_10mm_1925.mat

% Frex parameters 
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
Nfrex = length(frex);

freqbands = [frex(1) frex(9) frex(16) frex(21) frex(25) frex(30)];
freqnames = {'delta', 'theta', 'alpha', 'lowbeta', 'highbeta'};

Nvoxin = length(find(source.inverse.inside));

%% Avgs
load([p.results '\delta\voxshape_s1.mat']) % Load a sample to obtain filenames
features = fieldnames(voxshape{1});

trls = zeros(length(freqnames), Nsub, Nvoxin); 
avgshape = NaN(Nsub, Nvoxin,length(features));
for fb = 1:length(freqnames)
    frqband = freqnames{fb};
    display(frqband)
    for s = 1:Nsub
        load([p.results '\' frqband '\voxshape_s' num2str(s) '.mat'])
        for vx = 1:length(voxshape)
            trls(fb,s,vx) = length(voxshape{vx}.band_amp);
                voxshape{vx}.time_rdsym = abs(voxshape{vx}.time_rdsym-0.5)*2; % 0.5 [0 1]. Rescale to average [0 1]
                voxshape{vx}.time_ptsym = abs(voxshape{vx}.time_ptsym-0.5)*2; % 0.5 [0 1]. Rescale to average [0 1]
                avgshape(s,vx,:) = structfun(@(x) median(double(x)), voxshape{vx});             
        end
    end
    filename = p.results + "\avgshapes_" + freqnames{fb} + ".mat";
    save(filename, 'avgshape');
end


%% Statistics: ttest
addpath('Z:\Toolbox\slanCM\slanCM')
addpath('Z:\Toolbox\fieldtrip-20230118') 
freqnames = {'delta', 'theta', 'alpha', 'lowbeta', 'highbeta'};
paramlist = {'period', 'time_rdsym', 'time_ptsym', 'band_amp'};
[~, paramidx] = ismember(paramlist, features);
cmap = slanCM('vik', 256);
cmap(100:156, :) = [];
cmap_no_white = interp1(linspace(0, 1, size(cmap,1)), cmap, linspace(0, 1, 256));
mincolaxis = zeros(length(freqnames), length(paramlist));
maxcolaxis = zeros(length(freqnames), length(paramlist));

for fb = 1:length(freqnames)
    mkdir([p.results '/' freqnames{fb}])
    frqband = freqnames{fb};
    load([p.results '\avgshapes_' freqnames{fb}])
    for param = 1:length(paramlist)
        % Load avg data
        tmpdata = squeeze(avgshape(:,:,paramidx(param)));
        zData = normalize(tmpdata, 2, 'zscore');

        % Ttest
        [~, pvals,~, stats] = arrayfun(@(v) ttest(zData(:,v), 0), 1:size(tmpdata, 2),'UniformOutput',false);
        pvals = cell2mat(pvals);
        tVals = cellfun(@(s) s.tstat, stats);
        pfdr = mafdr(pvals(:), 'BHFDR', true);
        tValsfdr = tVals;
        tValsfdr(pfdr>.05) = 0;
        
        % Permutation Based Analysis
        nperm = 1000;
        maxTval = zeros(nperm, 1);
        for pp = 1:nperm
            permsign = (rand(Nsub, 1) > 0.5)*2 - 1; 
            permdata = zData .* permsign;
            [~,~,~, permstats] = arrayfun(@(v) ttest(permdata(:,v), 0), 1:size(permdata, 2),'UniformOutput',false);
            permtVals = cellfun(@(s) s.tstat, permstats);
            maxTval(pp) = max(abs(permtVals));
        end

        permTcorr = arrayfun(@(t) (sum(maxTval >= t) + 1) / (nperm + 1), abs(tVals));
        tValspermcor = tVals;
        tValspermcor(find(permTcorr > 0.05)) = 0;

        % sBOSC_nii(tValspermcor,[p.figures '/' freqnames{fb} '/Perm12_' paramlist{param} '_' frqband])

        maxv = max(tVals);
        minv = min(tVals);
        minlim = min(abs([maxv minv]));
        cfg = [];
        cfg.colmap = cmap_no_white;
        cfg.colim = [-10 10];
        cfg.interp = 'linear';
        sBOSC_sourcefig(tValspermcor, cfg)

        minval = min(tValspermcor(tValspermcor>0));
        if ~isempty(minval)
            mincolaxis(fb,param) = minval;
        end

        maxval = max(tValspermcor(tValspermcor<0));
        if ~isempty(maxval)
            maxcolaxis(fb,param) =  maxval;
        end
       
        close all
    end
end

%% ROIs
load(['Z:\OMEGA\OMEGA-NaturalFrequencies-main\mat_files\aal_voxel_label_10mm.mat'])
Nroi = length(aal_label_reduc);
Nbands = length(freqnames);

[~, idx_rdsym] = ismember('time_rdsym', features);
[~, idx_ptsym] = ismember('time_ptsym', features);

rdsym_roi = nan(Nsub, Nroi, Nbands);
ptsym_roi = nan(Nsub, Nroi, Nbands);
for fb = 1:Nbands
    load([p.results '\avgshapes_' freqnames{fb}])  
    
    rdsym_mat = squeeze(avgshape(:, :, idx_rdsym));  
    ptsym_mat = squeeze(avgshape(:, :, idx_ptsym)); 
    
    for roi = 1:Nroi
        inds = find(label_inside_aal_reduc == roi);
        voxs = voxel_inside_aal(inds);
        rdsym_roi(:, roi, fb) = mean(rdsym_mat(:, voxs), 2, 'omitnan');
        ptsym_roi(:, roi, fb) = mean(ptsym_mat(:, voxs), 2, 'omitnan');
    end
end

   
Ncols = 8;
Nrows = 5;
piecol = jet_omega_mod;
piecols = [piecol(1,:); piecol(20,:); piecol(40,:); piecol(55,:); piecol(64,:)];
aal_label_reduc = {'PreCG', 'SFG', 'SFG (orb)', 'MFG', 'MFG (orb)', 'IFG (oper)', 'IFG (tri)', 'IFG (orb)', 'Rolandic', 'SMA', 'Olfactory', 'SFG (med)', 'SFG (med orb)', 'Rectus', 'Insula', 'ACG', 'MCG', 'PCG', 'Hippoc', 'ParaHippoc', 'Calcarine', 'Cuneus', 'Lingual', 'SOG', 'MOG', 'IOG', 'Fusiform', 'PostCG', 'SPG', 'IPG', 'Supramarg', 'Angular', 'Precuneus', 'ParaCL', 'Heschl', 'STG', 'Temp pole (sup)', 'MTG', 'Temp pole (mid)', 'ITG'};
aal_label_sort = {'PreCG', 'SFG', 'SFG (orb)', 'MFG', 'MFG (orb)', 'IFG (oper)', 'IFG (orb)',  'IFG (tri)', 'Rolandic', 'PostCG', 'SPG', 'IPG', 'Supramarg', 'Angular', 'PCG', 'ParaCL', 'Hippoc', 'ParaHippoc', 'Temp pole (mid)', 'Temp pole (sup)', 'STG', 'MTG', 'Heschl', 'ITG',  'Precuneus', 'Calcarine', 'Cuneus', 'Lingual', 'SOG', 'MOG', 'IOG', 'Fusiform', 'SMA', 'Olfactory', 'SFG (med)', 'SFG (med orb)', 'Rectus', 'Insula', 'ACG', 'MCG'};
[~, sorted_aal] = ismember(aal_label_sort, aal_label_reduc);
aal_label_reduc = aal_label_reduc(sorted_aal);
rdsym_roi       = rdsym_roi(:,sorted_aal,:);
ptsym_roi       = ptsym_roi(:,sorted_aal,:);
for fb = 1:Nbands
    figure('Position', [100 100 Ncols*180 Nrows*180], 'Color', 'w');
    tl = tiledlayout(Nrows, Ncols, 'TileSpacing', 'compact', 'Padding', 'tight');

    for roi = 1:Nroi
        x = rdsym_roi(:, roi, fb);
        valid = ~isnan(x);

        nexttile;
        histogram(x(valid), 10, 'FaceColor', piecols(fb,:), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.7, 'Normalization', 'probability');
        hold on;
        xline(median(x(valid), 'omitnan'), 'k-', 'LineWidth', 1.5);
        xline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);  % 0.5 = symmetry
        title(aal_label_reduc{roi}, 'Interpreter', 'none', 'FontSize', 10);
        set(gca, 'FontSize', 2.5, 'Box', 'off', 'TickDir', 'out');
        xlim([0.15 0.35])
        hold off;
    end

        fig = gcf;
    fig.Units = 'centimeters';
    fig.Position(3) = 9;    % width = 9 cm
    fig.Position(4) = 5;    % adjust height proportionally

    exportgraphics(fig, fullfile(p.figures, ['rdsym_dist_roi_' freqnames{fb} '.png']), ...
        'Resolution', 300);

    close;

    % Same for ptsym
    figure('Position', [100 100 Ncols*180 Nrows*180], 'Color', 'w');
    tl = tiledlayout(Nrows, Ncols, 'TileSpacing', 'compact', 'Padding', 'tight');

    for roi = 1:Nroi
        x = ptsym_roi(:, roi, fb);
        valid = ~isnan(x);

        nexttile;
        histogram(x(valid), 10, 'FaceColor', piecols(fb,:), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.7, 'Normalization', 'probability');
        hold on;
        xline(median(x(valid), 'omitnan'), 'k-', 'LineWidth', 1.5);
        xline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        title(aal_label_reduc{roi}, 'Interpreter', 'none', 'FontSize', 10);
        set(gca, 'FontSize', 2.5, 'Box', 'off', 'TickDir', 'out');
        xlim([0.15 0.25])
        hold off;
    end

    fig = gcf;
    fig.Units = 'centimeters';
    fig.Position(3) = 9;    % width = 9 cm
    fig.Position(4) = 5;    % adjust height proportionally

    exportgraphics(fig, fullfile(p.figures, ['ptsym_dist_roi_' freqnames{fb} '.png']), ...
        'Resolution', 300);

    close;
end


%% Raw values
raw_group_means = cell(length(freqnames), length(paramlist));

for fb = 1:length(freqnames)
    frqband = freqnames{fb};
    
    load([p.results '\avgshapes_' freqnames{fb}])
    
    for param = 1:length(paramlist)
        tmpdata = squeeze(avgshape(:,:,paramidx(param)));        
        mean_raw_data = mean(tmpdata, 1, 'omitnan');
        raw_group_means{fb, param} = mean_raw_data;
        
     end
end

mean_period_alpha = raw_group_means{3,1};
sBOSC_nii(mean_period_alpha, [p.figures '\mean_period'])

%% Raw values to NIfTI (Frequency in Hz)
fsample = 256; 
mean_period_samples = raw_group_means{3, 1};
mean_freq_alpha_hz = fsample ./ mean_period_samples;
sBOSC_nii(mean_freq_alpha_hz, fullfile(p.figures, 'mean_freq_alpha_hz'));