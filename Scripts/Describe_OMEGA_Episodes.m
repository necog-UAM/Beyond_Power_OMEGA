%% sBOSC Oscillatory Episode Analysis
% Loads source-space oscillatory episodes detected with sBOSC (>=3 cycles)
% and characterizes their prevalence, duration, and power across frequency
% bands and brain regions.


%% Paths
p.data = 'Z:\OMEGA\OMEGA_data';
p.results = 'Z:\OMEGA\Enrique\Results';
p.figures = 'Z:\OMEGA\Enrique\figures';

% Toolbox paths
addpath('Z:\Toolbox\fieldtrip-20230118') % Path to Fieldtrip
ft_defaults
addpath(genpath('Z:\Toolbox\sourceBOSC'))
addpath('Z:\Toolbox\slanCM\slanCM') % colors to plot

srate = 128;

%% OMEGA 
% Participant and session list
sub = [1 2 3 4 5 6 7 8 9 11 12 14 15 16 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 37 39 40 41 42 44 45 46 47 48 49 50 51 52 55 56 57 58 59 60 61 62 63 64 65 67 68 69 70 71 72 73 74 75 76 77 78 79 80 84 85 87 88 89 90 91 92 94 95 96 97 98 99 101 102 103 104 105 106 134 145 146 148 149 150 151 152 154 155 156 157 158 159 160 161 165 166 167 168 169 170 171 175 176 177 179 181 184 185 195 197 200 207 208 210 212]';
ses = [1 1 1 1 1 1 1 2 1  2  1  2  1  1  1  2  3  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  3  1  2  1  2  1  2  1  1  2  1  1  1  1  1  1  1  1  2  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 ]';

subs = arrayfun(@(x) sprintf('%04d', x), sub, 'UniformOutput', false);
sess = arrayfun(@(x) sprintf('%04d', x), ses, 'UniformOutput', false);

nSub  = length(sub);

% Log-spaced frequencies from ~1.8 Hz to ~40 Hz
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
nFrex = length(frex);

freqbands = [frex(1) frex(9) frex(16) frex(21) frex(25) frex(30)];
freqnames = {'delta', 'theta', 'alpha', 'lowbeta', 'highbeta'};
nBands = length(freqnames);

load source_template_10mm_1925.mat % In sBOSC. Template with the voxels used (1925)
nVoxin = length(find(source.inverse.inside));

%% Preallocate
% Global % of recording time with at least voxel oscillating at any frequency
osc_global    = zeros(nSub, 1);
% Per voxel, collapsed across all frequencies
osc_any_freq  = zeros(nSub, nVoxin);
% Per voxel × frequency (up to 33 Hz)
osc_by_freq   = zeros(nSub, nVoxin, nFrex - 2);
% Per voxel × frequency band
osc_by_band   = zeros(nSub, nVoxin, nBands);

% Temporal overlap between band pairs
band_overlap   = zeros(nSub, nVoxin, nBands, nBands);
% Exclusive oscillation, % time oscillating in one band only
exclusive_band = zeros(nSub, nVoxin, nBands);
 
% Episode properties 
all_dur_secs_by_band = cell(nSub, nVoxin, nBands); % Duration in seconds
all_dur_cyc_by_band  = cell(nSub, nVoxin, nBands); % Duration in cycles
all_pow_by_band      = cell(nSub, nVoxin, nBands); % Mean power per episode

%%  Part I. Extract episode properties
for s=1:nSub
    s
    load([p.data '\sub-' subs{s} '\ses-' sess{s} '\conepis\epis_v1.mat']) % to check size of time points
    nTp = length(datasource.time);
    load([p.data '\sub-' subs{s} '\ses-' sess{s} '\episorig3cyc.mat']) % Variable with only epis > 3 cycle


%% Transform epis from cells to matrix and extract power and duration
    episodes = false(nVoxin,nFrex,nTp);

    temp_dur_sec = cell(nVoxin, length(freqnames)); 
    temp_dur_cyc = cell(nVoxin, length(freqnames)); 
    temp_pow     = cell(nVoxin, length(freqnames));
    
    for vx = 1:nVoxin % Voxel loop
        for ep = 1:length(epis3c{vx}) % Episode loop
            fm = epis3c{vx}(ep).freq;
            fbin = dsearchn(frex',fm);
            tpts = epis3c{vx}(ep).timeps;

            episodes(vx,fbin(1),tpts) = true;
            
            % EPisode duration
            dur_sec = length(tpts) / srate; % in seconds
            dur_cyc = dur_sec * fm; % in cycles

            % Episode power
            mean_pow = mean(epis3c{vx}(ep).power, 'omitmissing');

            % Extract freq band of this episode
            fb_idx = find(fm >= freqbands(1:end-1) & fm <= freqbands(2:end));
            
            if ~isempty(fb_idx)
                fb_target = fb_idx(1);
                temp_dur_sec{vx, fb_target} = [temp_dur_sec{vx, fb_target}, dur_sec];
                temp_dur_cyc{vx, fb_target} = [temp_dur_cyc{vx, fb_target}, dur_cyc];  
                temp_pow{vx, fb_target}     = [temp_pow{vx, fb_target}, mean_pow]; 
            end
        end
    end
    
    % Cut at 33 Hz
    episodes = episodes(:,1:30,:);
    
    % Store
    all_dur_secs_by_band(s, :, :) = temp_dur_sec;
    all_dur_cyc_by_band(s, :, :)  = temp_dur_cyc;
    all_pow_by_band(s, :, :)     = temp_pow;


%% Oscillatory prevalence
        
        % Percent of  where at least one voxel oscillates at any freq
        any_tp = squeeze(sum(sum(episodes, 2), 1) > 0);
        osc_global(s) = sum(any_tp) ./ nTp * 100;

        % Percent of time points with any frequency per voxel
        any_freq_vox = squeeze(sum(episodes, 2) > 0);
        osc_any_freq(s,:) = sum(any_freq_vox, 2) ./ nTp * 100;

        % Percent per voxel and per frequency  
        osc_by_freq(s,:,:) = sum(episodes, 3) ./ nTp * 100;

        % Percent per frequency band
        fb_mask = false(nVoxin, nBands, nTp);

        for fb = 1:nBands
            if fb == nBands
                ct = 0; %
            else
                ct = 1;
            end
            fband = [freqbands(fb) freqbands(fb+1)];
            fbidx = dsearchn(frex', fband');
            fb_mask(:, fb, :) = squeeze(sum(episodes(:, fbidx(1):(fbidx(2)-ct), :), 2)) > 0;
        end

        % Band overlaps
        overlap_voxels = zeros(nVoxin, nBands, nBands);

        for fb = 1:nBands
            thisband = squeeze(fb_mask(:, fb, :));
            
            % Prevalence within this band
            osc_by_band(s, :, fb) = sum(thisband, 2) ./ nTp * 100;
            
            % Overlaps
            for fb2 = 1:nBands
                freq_overlap = squeeze(fb_mask(:, fb2, :));
                overlap = thisband .* freq_overlap; 
                overlap_voxels(:, fb, fb2) = sum(overlap, 2) ./ nTp * 100;
            end
            
            % Exclusive. Active in this band but not in others
            other_idx   = setdiff(1:Nbands, fb);
            other_bands = squeeze(sum(fb_mask(:, other_idx, :), 2)) > 0;
            exclusive_band(s, :, fb) = sum(thisband & ~other_bands, 2)' / nTp * 100;        
        end
        
        band_overlap(s, :, :, :) = overlap_voxels;
                    
end

% Store results
save(fullfile(p.results, 'oscillatory_results_3cyc.mat'), 'osc_global', 'osc_any_freq', 'osc_by_band','osc_by_freq','band_overlap','exclusive_band', 'all_dur_secs_by_band', 'all_dur_cyc_by_band', 'all_pow_by_band');


%% Part II. Analysis and figures
load(fullfile(p.results, 'oscillatory_results_3cyc.mat'))

% Global oscillatory prevalence
    mean_osc_global = mean(osc_global);
    sd_osc_global = std(osc_global, 1);

% Average of oscillatory prevalence per voxel across frequencies
    mean_osc_any_freq = squeeze(mean(osc_any_freq, 1)); 
    sd_osc_any_freq = squeeze(std(osc_any_freq, 0, 1));

        % Figure and .nii
        cfg = [];
        cfg.colmap = slanCM('BuPu');
        cfg.interp = 'linear';
        cfg.colim = [0 60];
        sBOSC_sourcefig(mean_osc_any_freq, cfg)
        sBOSC_nii(mean_osc_any_freq, [p.figures '\mean_osc_any_freq'])
     
% Prevalence in frequnecy bands
    colmap = {'Blues', 'Greens', 'YlOrBr', 'Oranges', 'Reds'};
    mean_osc_band = zeros(nVoxin, nBands);

        % Figure and .nii
        for fb = 1:nBands
                mean_osc_band(:, fb) = squeeze(mean(osc_by_band(:, :, fb), 1));       
                cfg = [];
                cfg.colmap = slanCM(colmap{fb});
                cfg.savefig = 'no';
                cfg.interp = 'spline';
                sBOSC_sourcefig(mean_osc_band(:, fb), cfg); 
                sgtitle(['avg_' freqnames{fb} '_osc']);
                sBOSC_nii(mean_osc_band(:, fb), ['avg_' freqnames{fb} '_osc']);
        end


% Global oscillatory and non-oscillatory time
    osc = mean(osc_any_freq, 'all');
    nonosc = 100 - osc;

    % Figure
        figure, pie([nonosc osc], {'Non-oscillatory', 'Oscillatory'})
        colormap("hot")
    
    mean_per_subject = mean(osc_any_freq, 2); 
    sd_across_participants = std(mean_per_subject);
    mean_per_voxel = mean(osc_any_freq, 1);
    sd_across_voxels = std(mean_per_voxel);
 
% Frequency bands
    band_global_means = squeeze(mean(mean(osc_by_band, 1), 2)); 
    band_sectors = (band_global_means ./ osc) .* 100;

% Frequency band overlap and exclusive
    mean_band_overlap = squeeze(mean(band_overlap, [1, 2]));
    mean_exclusive_band = squeeze(mean(exclusive_band, [1, 2]));

% Normalize to the global oscillatory time
    overlap_graph = (mean_band_overlap / osc) * 100;
    exclusive_graph = (mean_exclusive_band / osc) * 100;
    exclusive_graph ./  diag(overlap_graph) * 100 % Percent of time a band is oscillating with no concurrent bands

        % Figure
        piecol = jet_omega_mod; % colormap from sBOSC
        piecols = [piecol(1,:); piecol(20,:); piecol(40,:); piecol(55,:); piecol(64,:)];
    
        X_coords = [0, 2, 0, -2, 0, 1.5];
        Y_coords = [0, 0, 2, 0, -2, 1.5];
    
        scaling_factor = 500;
        exclusive_color = [0.7 0.7 0.7]; 
            
        idx_fb = 1:nBands;
        
        for fb = 1:nBands
                figure;
                percentages = zeros(1, nBands + 1);
                
                tmp = overlap_graph(fb, :);
                percentages(1) = tmp(fb); % Self-overlap
                
                idx2_fb = idx_fb(idx_fb ~= fb); % All other bands
                percentages(2:nBands) = tmp(idx2_fb); % Overlap with other bands
                percentages(nBands + 1) = exclusive_graph(fb); 
                
                sizes = percentages * scaling_factor;     
                Colors = [piecols(fb, :); piecols(idx2_fb, :); exclusive_color];
                
                scatter(X_coords, Y_coords, sizes, Colors, 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.6);
             
                for i = 1:length(percentages)
                    pp = percentages(i);
                    x_pos = X_coords(i);
                    y_pos = Y_coords(i);
                    
                    if pp < 50
                        text_color = [0 0 0]; 
                    else
                        text_color = [1 1 1]; 
                    end
                    
                    if i == (nBands + 1)
                        label_str = sprintf('%.2f%%\n(Only)', pp);
                    else
                        label_str = sprintf('%.2f%%', pp);
                    end
                    
                    text(x_pos, y_pos, label_str, ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'FontSize', 14, ...
                        'Color', text_color, ...
                        'FontWeight', 'bold');
                end
                
                axis([-2.5 2.5 -2.5 2.5])
                axis equal;
                axis off;
                
                filename = fullfile(p.figures, [freqnames{fb} '_overlap']);
                saveas(gcf, filename, 'png');
                close;
            end


%% Part III ROI Analysis

    load(['Z:\Enrique\Descriptives\Files\aal_voxel_label_10mm.mat']) %%% in Files
    Nroi = length(aal_label_reduc);
    
    % Median prevalence per ROI and band
    rois = zeros(Nroi, length(freqnames), nSub);
    roivoxs = zeros(nVoxin, length(freqnames), nSub);
    for roi = 1:Nroi
        for fb = 1:length(freqnames)
        inds = find(label_inside_aal_reduc==roi);
        voxs = voxel_inside_aal(inds);
        rois(roi,fb,:) = squeeze(median(osc_by_band(:,voxs,fb),2));
        roivoxs(voxs,fb,:) = repmat(squeeze(median(osc_by_band(:,voxs,fb),2))',length(voxs),1);
        end
    end
    
        % Figure
        Ncols = 8;
        Nrows = 5;
        
        figure
        t = tiledlayout(Nrows, Ncols, 'TileSpacing', 'compact','Padding', 'tight');
        for roi = 1:Nroi
            boxdata = [squeeze(rois(roi,1,:)) squeeze(rois(roi,2,:)) squeeze(rois(roi,3,:)) squeeze(rois(roi,4,:)) squeeze(rois(roi,5,:))];
            nexttile(roi);
            boxplot(boxdata, 'Symbol', '', 'Widths', 0.8, 'Color', 'k', 'Notch', 'off');
            h = flipud(findobj(gca, 'Tag', 'Box'));
        
                for k = 1:length(h)
                    XData = get(h(k), 'XData');
                    YData = get(h(k), 'YData');    
                    patch(XData, YData, piecols(k,:), 'FaceAlpha', 0.7);    
                    delete(h(k));
                end
        
            Q1 = prctile(boxdata, 25);
            Q3 = prctile(boxdata, 75);
            IQR = Q3 - Q1;
            lower = Q1 - 1.5*IQR;
            lower(lower<0) = 0;
            upper = Q3 + 1.5*IQR;
            buffer = max(IQR)*0.2;
            ylim([min(lower)-buffer max(upper)+buffer]);
            
            title([aal_label_reduc{roi}], 'FontSize', 8, 'Interpreter', 'none'); 
            xticklabels([]);
            set(gca, 'FontSize', 7);
            set(gca, 'YGrid', 'on', 'GridAlpha', 0.4, 'GridColor', [0.6 0.6 0.6]);
            box off;
        set(gca, 'TickDir', 'out');
        set(gca, 'TickLength', [0.005 0.005]);
        end
        linkaxes(t.Children, 'y');

% ROIs prevalence across frequencies (prevalence per ROI at any freq)

    roivoxind = zeros(1,nVoxin);
    rois_any_freq = zeros(Nroi, nSub);
    for roi = 1:Nroi
        inds = find(label_inside_aal_reduc==roi);
        voxs = voxel_inside_aal(inds);
        rois_any_freq(roi, :) = squeeze(median(osc_any_freq(:,voxs),2));
        roivoxind(voxs) = roi;
    end

    roi_mean_any_freq = mean(rois_any_freq, 2);
    roi_std_any_freq  = std(rois_any_freq, 0, 2);
    
        % Table
        ROI_Table_AnyFreq = table(aal_label_reduc(:), roi_mean_any_freq, roi_std_any_freq, 'VariableNames', {'ROI_Name', 'AnyFreq_Mean_Pct', 'AnyFreq_Std_Pct'});
        writetable(ROI_Table_AnyFreq, fullfile(p.results, 'ROI_AnyFreq_Percentages.csv'));


% ROIs prevalence by frequency band

        % Table
        ROI_Table = table(aal_label_reduc(:), 'VariableNames', {'ROI_Name'});
        for fb = 1:Nbands
            raw_band_data = squeeze(rois(:, fb, :));
            ROI_Table.(sprintf('%s_Mean_Pct', freqnames{fb})) = mean(raw_band_data, 2);
            ROI_Table.(sprintf('%s_Std_Pct',  freqnames{fb})) = std(raw_band_data,  0, 2);
        end
        writetable(ROI_Table, fullfile(p.results, 'ROI_Raw_Percentages.csv'));


% ROI z-scored across voxels
    
    z_osc_by_band = zeros(nSub, nVoxin, nBands);
    
    for s = 1:nSub
        for fb = 1:nBands
            curr_data = osc_by_band(s, :, fb);
            z_osc_by_band(s, :, fb) = (curr_data - mean(curr_data)) ./ std(curr_data);
        end
    end
    
    rois_z = zeros(Nroi, nBands, nSub);
    
    for roi = 1:Nroi
        inds = find(label_inside_aal_reduc == roi);
        voxs = voxel_inside_aal(inds);
        for fb = 1:nBands
            rois_z(roi, fb, :) = squeeze(median(z_osc_by_band(:, voxs, fb), 2));
        end
    end

        % Figure
        figure('Name', 'ROI z-scores by Band', 'Color', 'w', 'Position', [50 50 1400 800]);
        t = tiledlayout(Nrows, Ncols, 'TileSpacing', 'compact', 'Padding', 'tight');     
        for roi = 1:Nroi
            boxdata = squeeze(rois_z(roi, :, :))'; 
            nexttile(roi);
            boxplot(boxdata, 'Symbol', '', 'Widths', 0.8, 'Color', 'k', 'Notch', 'off');
         
            h = flipud(findobj(gca, 'Tag', 'Box'));
            for k = 1:length(h)
                patch(get(h(k), 'XData'), get(h(k), 'YData'), piecols(k, :), 'FaceAlpha', 0.7);
                delete(h(k));
            end
         
            hold on;
            yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
         
            title(aal_label_reduc{roi}, 'FontSize', 8, 'Interpreter', 'none');
            xticklabels([]);
            ylim([-2.5 4]);
            set(gca, 'FontSize', 7, 'YGrid', 'on', 'GridAlpha', 0.3, 'Box', 'off');
            hold off;     
        end     
        linkaxes(t.Children, 'y');

% Brain maps: areas 

    idx_frontal   = 1:8; 
    idx_parietal  = [9,18, 28:32 34]; 
    idx_temporal  = [19:20, 35:40]; 
    idx_occipital = [21:27, 33]; 
    idx_medial    = 10:17;
    
    vox_frontal   = double(ismember(roivoxind, idx_frontal));
    vox_parietal  = double(ismember(roivoxind, idx_parietal));
    vox_temporal  = double(ismember(roivoxind, idx_temporal));
    vox_occipital = double(ismember(roivoxind, idx_occipital));
    vox_medial    = double(ismember(roivoxind, idx_medial));
    
    darkcol = [0.4, 0.0, 0.6]; 
    lightcol  = [0.7, 0.4, 0.8];
    ncolor = 64;
    r = linspace(lightcol(1), darkcol(1), ncolor)';
    g = linspace(lightcol(2), darkcol(2), ncolor)';
    b = linspace(lightcol(3), darkcol(3), ncolor)';
    purplemap = [r, g, b];
    
    area_masks = struct( ...
        'vox_frontal',   vox_frontal,   ...
        'vox_parietal',  vox_parietal,  ...
        'vox_temporal',  vox_temporal,  ...
        'vox_occipital', vox_occipital, ...
        'vox_medial',    vox_medial);
    area_names = fieldnames(area_masks);

        % Figure
        cfg        = [];
        cfg.colmap = purplemap;
        cfg.colim  = [-1 1];
        cfg.interp = 'spline';

        for ar = 1:length(area_names)
            name = area_names{ar};
            sBOSC_sourcefig(area_masks.(name), cfg);        
        end


%% Part IV ANOVAs


area_idx   = {find(vox_frontal), find(vox_parietal), find(vox_temporal), find(vox_occipital), find(vox_medial)};
area_names = {'Frontal', 'Parietal', 'Temporal', 'Occipital', 'Medial'};
nAreas     = length(area_names);
area_data = nan(nSub, nBands, nAreas);
for ar = 1:nAreas
    voxs = area_idx{ar};
    for fb = 1:nBands
        area_data(:, fb, ar) = mean(z_osc_by_band(:, voxs, fb), 2, 'omitnan');
    end
end

%% Repeated-measures ANOVA. Prevalence band diffs within each area

    band_pairs  = nchoosek(1:nBands, 2);
    rm_formula = sprintf('%s-%s~1', freqnames{1}, freqnames{end});
    ranova_results = cell(nAreas, 1);
    mc_results     = cell(nAreas, 1);
  
        for ar = 1:nAreas
            region_data = squeeze(area_data(:, :, ar));
            T_rm        = array2table(region_data, 'VariableNames', freqnames);
            within      = table(freqnames', 'VariableNames', {'Band'});
            rm             = fitrm(T_rm, rm_formula, 'WithinDesign', within);
            ranova_out     = ranova(rm, 'WithinModel', 'Band');
            mc             = multcompare(rm, 'Band', 'ComparisonType', 'bonferroni');
            ranova_results{ar} = ranova_out;
            mc_results{ar}     = mc;
        end



        % Figure
        figure('Color', 'w', 'Position', [100 100 1200 350]);
        tl = tiledlayout(1, Nregions, 'TileSpacing', 'compact', 'Padding', 'tight'); 
        letters  = 'abcde';
        y_letter = 1.5;
         
        for ar = 1:nAreas
            ardata = squeeze(area_data(:, :, ar));
            mc          = mc_results{ar};
            letter_label = letters(1:Nbands);
            for pair = 1:size(band_pairs, 1)
                b1  = freqnames{band_pairs(pair, 1)};
                b2  = freqnames{band_pairs(pair, 2)};
                row = (strcmp(mc.Band_1, b1) & strcmp(mc.Band_2, b2)) | (strcmp(mc.Band_1, b2) & strcmp(mc.Band_2, b1));
                if mc.pValue(row) >= 0.05
                    letter_label(band_pairs(pair, 2)) = letter_label(band_pairs(pair, 1));
                end
            end   
            nexttile;
            boxplot(ardata, freqnames, 'Symbol', '', 'Colors', 'k');
            h = flipud(findobj(gca, 'Tag', 'Box'));
            for k = 1:length(h)
                patch(get(h(k), 'XData'), get(h(k), 'YData'), piecols(k, :), 'FaceAlpha', 0.7);
                delete(h(k));
            end
            hold on;
            yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
            for k = 1:Nbands
                text(k, y_letter, letter_label(k), ...
                    'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
            end
            ylim([-1.75 1.75]);
            title(area_names{ar}, 'FontSize', 11);
            if ar == 1, ylabel('Z-score', 'FontSize', 11); end
            set(gca, 'FontSize', 10, 'Box', 'off', 'YGrid', 'on', ...
                'GridAlpha', 0.3, 'TickDir', 'out');
            hold off;
        end 
        exportgraphics(gcf, fullfile(p.figures, 'boxplot_anova_rm_areas.png'), 'Resolution', 300);


%% Episode duration

    mean_dur_cyc_subj = nan(nSub, Nbands);
    mean_dur_sec_subj = nan(nSub, Nbands);
     
    for s = 1:nSub
        for fb = 1:Nbands
     
            cyc_cells = all_dur_cyc_by_band(s, :, fb);
            sec_cells = all_dur_secs_by_band(s, :, fb);
    
            cyc_cells = cyc_cells(~cellfun(@isempty, cyc_cells));
            sec_cells = sec_cells(~cellfun(@isempty, sec_cells));
     
            if ~isempty(cyc_cells)
                all_cyc = cell2mat(cellfun(@(x) x(:), cyc_cells(:), 'UniformOutput', false));
                all_sec = cell2mat(cellfun(@(x) x(:), sec_cells(:), 'UniformOutput', false));
     
                mean_dur_cyc_subj(s, fb) = mean(all_cyc);
                mean_dur_sec_subj(s, fb) = mean(all_sec);
            end
     
        end
    end
     
    % Grand means across participants
    grand_dur_cyc_mean = mean(mean_dur_cyc_subj, 1, 'omitnan');
    grand_dur_cyc_sd   = std( mean_dur_cyc_subj, 0, 1, 'omitnan');
    grand_dur_sec_mean = mean(mean_dur_sec_subj, 1, 'omitnan');
    grand_dur_sec_sd   = std( mean_dur_sec_subj, 0, 1, 'omitnan');


%% Repeated-measures ANOVA.Duration band diffs

    T_rm       = array2table(mean_dur_cyc_subj, 'VariableNames', freqnames);
    within     = table(freqnames', 'VariableNames', {'Band'});
    rm_formula = sprintf('%s-%s~1', freqnames{1}, freqnames{end});
    rm         = fitrm(T_rm, rm_formula, 'WithinDesign', within);
    ranova_out = ranova(rm, 'WithinModel', 'Band');
    mc         = multcompare(rm, 'Band', 'ComparisonType', 'bonferroni');
    
    % Post-hoc t-statistics and p-values
    t_vals = mc.Difference ./ mc.StdErr;
    p_vals = mc.pValue;
     
    results_table = table(mc.Band_1, mc.Band_2, t_vals, p_vals, ...
        'VariableNames', {'Band1', 'Band2', 'tStat', 'pValue'});
    disp(results_table);

        % Figure
        figure('Color', 'w', 'Position', [100 100 500 400]);
        boxplot(mean_dur_cyc_subj, 'Labels', freqnames, 'Colors', 'k', 'Symbol', '.');
         
        h = flipud(findobj(gca, 'Tag', 'Box'));
        for k = 1:length(h)
            patch(get(h(k), 'XData'), get(h(k), 'YData'), piecols(k, :), 'FaceAlpha', 0.7);
            delete(h(k));
        end   
        ylabel('Duration (cycles)', 'FontSize', 12);
        set(gca, 'FontSize', 12, 'YGrid', 'on', 'GridAlpha', 0.3, 'Box', 'off', 'TickDir', 'out');

        % Histograms: cycles and seconds
        hist_params = struct( ...
            'data',     {all_dur_cyc_by_band,  all_dur_secs_by_band}, ...
            'binwidth', {1,                    0.2                  }, ...
            'xlim',     {[3 10],               [0 2]                }, ...
            'ylim',     {[0 0.6],              [0 0.8]              }, ...
            'xlabel',   {'Cycles',             'Seconds'            }, ...
            'name',     {'Duration — Cycles',  'Duration — Seconds' });
        
        for hp = 1:2
            figure('Name', ['Durations — ' hist_params(hp).name], ...
                'Position', [100 100 1200 300], 'Color', 'w');
            tiledlayout(1, Nbands, 'TileSpacing', 'compact');
            for fb = 1:Nbands
                band_cells  = hist_params(hp).data(:, :, fb);
                valid_cells = band_cells(~cellfun(@isempty, band_cells));
                all_bursts  = cell2mat(cellfun(@(x) x(:), valid_cells(:), 'UniformOutput', false));
                nexttile;
                histogram(all_bursts, 'BinWidth', hist_params(hp).binwidth, 'Normalization', 'probability', 'FaceColor', piecols(fb,:), 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                title(freqnames{fb}, 'Interpreter', 'none', 'FontSize', 14);
                xlim(hist_params(hp).xlim);
                ylim(hist_params(hp).ylim);
                xlabel(hist_params(hp).xlabel);
                if fb == 1, ylabel('Probability'); end
                set(gca, 'FontSize', 12, 'Box', 'off', 'TickDir', 'out');
            end        
        end


%% Part V Permutation Analysis. 

%% Duration and Power

dur_subj_vox = nan(nSub, Nvoxin, Nbands);
pow_subj_vox = nan(nSub, Nvoxin, Nbands);
 
for s = 1:nSub
    for vx = 1:Nvoxin
        for fb = 1:Nbands
            cyc_vals = all_dur_cyc_by_band{s, vx, fb};
            pow_vals = all_pow_by_band{s, vx, fb};
            if ~isempty(cyc_vals)
                dur_subj_vox(s, vx, fb) = mean(cyc_vals, 'omitnan');
                pow_subj_vox(s, vx, fb) = mean(pow_vals, 'omitnan');
            end
        end
    end
end

% Zscore across voxels
    dur_z_all = cell(Nbands, 1);
    pow_z_all = cell(Nbands, 1);
     
    for fb = 1:Nbands
        dur_mat = squeeze(dur_subj_vox(:, :, fb));   % [Nsub × Nvoxin]
        pow_mat = squeeze(pow_subj_vox(:, :, fb));
        dur_z_all{fb} = (dur_mat - mean(dur_mat, 2, 'omitnan')) ./ std(dur_mat, 0, 2, 'omitnan');
        pow_z_all{fb} = (pow_mat - mean(pow_mat, 2, 'omitnan')) ./ std(pow_mat, 0, 2, 'omitnan');
    end


cmap = slanCM('vik', 256);
cmap(100:156, :) = [];
cmap = interp1(linspace(0, 1, size(cmap,1)), cmap, linspace(0, 1, 256));
param_labels = {'dur_cyc', 'power'};
data_ep      = {dur_subj_vox, pow_subj_vox};
nperm        = 1000;
tstat_colmin = zeros(Nbands, length(param_labels));
tstat_colmax = zeros(Nbands, length(param_labels));

    for fb = 1:Nbands
        for param = 1:length(param_labels)
            tmpdata = squeeze(data_ep{param}(:, :, fb));
            zData   = (tmpdata - mean(tmpdata, 2, 'omitnan')) ./ std(tmpdata, 0, 2, 'omitnan');
            % observed t-stat against 0
            [~, ~, ~, stats] = arrayfun(@(v) ttest(zData(:, v), 0), 1:Nvoxin, 'UniformOutput', false);
            tVals = cellfun(@(s) s.tstat, stats);
     
            % Max-statistic permutation distribution
            maxT_perm = zeros(nperm, 1);
            for pp = 1:nperm
                permsign        = (rand(nSub, 1) > 0.5) * 2 - 1;
                [~,~,~, pstats] = arrayfun(@(v) ttest(zData(:, v) .* permsign, 0), 1:Nvoxin, 'UniformOutput', false);
                maxT_perm(pp)   = max(abs(cellfun(@(s) s.tstat, pstats)));
            end   
            % Permutation-corrected
            permP        = arrayfun(@(t) (sum(maxT_perm >= t) + 1) / (nperm + 1), abs(tVals));
            tVals_thresh = tVals;
            tVals_thresh(permP > 0.05) = 0;
             
            pos_vals = tVals_thresh(tVals_thresh > 0);
            neg_vals = tVals_thresh(tVals_thresh < 0);
            if ~isempty(pos_vals), tstat_colmin(fb, param) = min(pos_vals); end
            if ~isempty(neg_vals), tstat_colmax(fb, param) = max(neg_vals); end
     
            % Brain map
            cfg         = [];
            cfg.colmap  = cmap;
            cfg.colim   = [-10 10];
            cfg.interp  = 'linear';
            cfg.savefig = 'no';
            sBOSC_sourcefig(tVals_thresh, cfg);
            close all;
        end
    end


%% Correlation duration and power
    cmap_corr = flipud(slanCM('PuOr'));
    cmap_corr(100:156, :) = [];
    cmap_corr = interp1(linspace(0, 1, size(cmap_corr, 1)), cmap_corr, linspace(0, 1, 256));
    
    r_vox  = nan(nVoxin, nBands);
    t_vox  = nan(nVoxin, nBands);
    sig_vox = false(nVoxin, nBands);
    
    r_colmin = zeros(Nbands, 1);
    r_colmax = zeros(Nbands, 1);
    
    for fb = 1:nBands
        dur_z = dur_z_all{fb};
        pow_z = pow_z_all{fb};
        % Exclude voxels with insufficient valid participants
        n_valid   = sum(~isnan(dur_z) & ~isnan(pow_z), 1);
        valid_vox = n_valid >= 10;
        dur_c = dur_z(:, valid_vox);
        pow_c = pow_z(:, valid_vox);
        dur_c(isnan(dur_c)) = 0;
        pow_c(isnan(pow_c)) = 0;
        
        % Pearson r 
        num   = sum(dur_c .* pow_c, 1);
        denom = sqrt(sum(dur_c.^2, 1) .* sum(pow_c.^2, 1));
        r_vox(valid_vox, fb) = (num ./ denom)';
    
        % Convert r to t
        n = n_valid(valid_vox)';
        t_vox(valid_vox, fb) = r_vox(valid_vox, fb) .* sqrt((n - 2) ./ (1 - r_vox(valid_vox, fb).^2));
    
        % Max-statistic permutation 
        maxT_perm = zeros(nperm, 1);
        for pp = 1:nperm
            permsign    = (rand(nSub, 1) > 0.5) * 2 - 1;
            dur_perm    = dur_c .* permsign;
            num_p       = sum(dur_perm .* pow_c, 1);
            denom_p     = sqrt(sum(dur_perm.^2, 1) .* sum(pow_c.^2, 1));
            r_perm      = (num_p ./ denom_p)';
            t_perm      = r_perm .* sqrt((n - 2) ./ (1 - r_perm.^2));
            maxT_perm(pp) = max(abs(t_perm), [], 'omitnan');
        end
    
        % Permutation corrected
        t_obs   = t_vox(valid_vox, fb);
        permP_r = arrayfun(@(t) (sum(maxT_perm >= t) + 1) / (nperm + 1), abs(t_obs));
        sig_vox(valid_vox, fb) = permP_r < 0.05;
     
        r_sig = r_vox(:, fb);
        r_sig(~sig_vox(:, fb)) = 0;
        pos_r = r_sig(r_sig > 0); if ~isempty(pos_r), r_colmin(fb) = min(pos_r); end
        neg_r = r_sig(r_sig < 0); if ~isempty(neg_r), r_colmax(fb) = max(neg_r); end
     
        % Brain map 
        r_masked = r_vox(:, fb);
        r_masked(~sig_vox(:, fb)) = 0;
        cfg         = [];
        cfg.colmap  = cmap_corr;
        cfg.colim   = [-0.8 0.8];
        cfg.interp  = 'linear';
        cfg.savefig = 'no';
        sBOSC_sourcefig(r_masked, cfg);
        close all;
    end


%% Correlation duration and power (ROI)
    r_roi = nan(Nroi, Nbands);
    p_roi = nan(Nroi, Nbands);
    
    for fb = 1:Nbands
        dur_z = dur_z_all{fb};
        pow_z = pow_z_all{fb};
        for roi = 1:Nroi
            inds  = find(label_inside_aal_reduc == roi);
            voxs  = voxel_inside_aal(inds);
            x     = mean(dur_z(:, voxs), 2, 'omitnan');
            y     = mean(pow_z(:, voxs), 2, 'omitnan');
            valid = ~isnan(x) & ~isnan(y);
            if sum(valid) > 3
                [R_mat, P_mat]  = corrcoef(x(valid), y(valid));
                r_roi(roi, fb)  = R_mat(1, 2);
                p_roi(roi, fb)  = P_mat(1, 2);
            end
        end
     
    end

        % Figure
        NCols = 5;
         
        for fb = 1:Nbands
            dur_z = dur_z_all{fb};
            pow_z = pow_z_all{fb};
            % Significant ROIs sorted by r
            sig_rois  = p_roi(:, fb) < 0.05;
            r_band    = r_roi(:, fb);
            r_band(~sig_rois) = NaN;
            [~, idx_sorted] = sort(abs(r_band), 'descend', 'MissingPlacement', 'last');
         
            roi_plot  = idx_sorted(1:num_valid);
            num_rows  = ceil(num_valid / NCols);
         
            fig = figure('Position', [100 100 300*NCols 250*num_rows], 'Color', 'w');
            tiledlayout(num_rows, NCols, 'TileSpacing', 'compact', 'Padding', 'compact');
         
            for ri = 1:num_valid
                roi   = roi_plot(ri);
                inds  = find(label_inside_aal_reduc == roi);
                voxs  = voxel_inside_aal(inds);
                x     = mean(dur_z(:, voxs), 2, 'omitnan');
                y     = mean(pow_z(:, voxs), 2, 'omitnan');
                valid = ~isnan(x) & ~isnan(y);
                pv = p_roi(roi, fb);
                if pv < 0.001,     pstr = 'p < 0.001';
                elseif pv < 0.01,  pstr = 'p < 0.01';
                else,              pstr = 'p < 0.05';
                end
                % Regression line
                p_fit  = polyfit(x(valid), y(valid), 1);
                x_line = linspace(min(x(valid)), max(x(valid)), 100);
                nexttile;
                scatter(x(valid), y(valid), 30, piecols(fb, :), 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
                hold on; plot(x_line, polyval(p_fit, x_line), 'k-', 'LineWidth', 1.5);
                text(0.05, 0.90, sprintf('r = %.2f\n%s', r_roi(roi, fb), pstr), 'Units', 'normalized', 'FontSize', 16, 'Color', 'k');      
                title(strrep(aal_label_reduc{roi}, '_', ' '), 'Interpreter', 'none', 'FontSize', 20, 'FontWeight', 'normal');
         
                if ri > (num_rows - 1) * NCols, xlabel('Duration (Z)', 'FontSize', 10); end
                if mod(ri, NCols) == 1,         ylabel('Power (Z)',    'FontSize', 10); end

                set(gca, 'FontSize', 14, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1);
                hold off;   
            end    
            exportgraphics(fig, fullfile(p.figures, ...
                sprintf('scatter_sig_rois_%s.png', freqnames{fb})), 'Resolution', 300);
            close(fig);
        end
