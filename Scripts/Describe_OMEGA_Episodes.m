%% This script loads previously identified episodes with sBOSC with a duration >= 3 cycles.
%% Analyzes: prevalence, duration and power of oscillatory episodes.

%% Paths

p.data = 'Z:\OMEGA\OMEGA_data';
p.results = 'Z:\OMEGA\Enrique\Results';
p.figures = 'Z:\OMEGA\Enrique\figures';

% Toolbox paths
addpath('Z:\Toolbox\fieldtrip-20230118') % Path to Fieldtrip
addpath('Z:\Toolbox\NECOG')
ft_defaults

addpath(genpath('Z:\Toolbox\sourceBOSC'))
addpath('Z:\Toolbox\slanCM\slanCM') % colors to plot

addpath('Z:\OMEGA\OMEGA-NaturalFrequencies-main\subfxs') % para jet omega mod
srate = 128;

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
Nbands = length(freqnames);

Nvoxin = length(find(source.inverse.inside));

%% Part I: Extract data from episodes
osc_global = zeros(Nsub,1); % Percent of total time where at least one voxel is oscillating at any frequency [Sub x 1]
osc_any_freq = zeros(Nsub,Nvoxin); % Percent of time each voxel is oscillating (collapsed across all frequencies). [Sub x Vox]
osc_by_freq = zeros(Nsub,Nvoxin, length(frex)-2); % Percent of time each voxel is oscillating at each frequency. [Sub x Vox x Freq]
osc_by_band = zeros(Nsub,Nvoxin, 5); % Percent of time each voxel is oscillating within frequency bands. [Sub x Vox x Band]
band_overlap =  zeros(Nsub,Nvoxin,length(freqnames),length(freqnames)); % Percent of temporal overlap in oscillations between pairs of bands per voxel. [Sub x Vox x Band1 x Band2]
exclusive_band = zeros(Nsub,Nvoxin, length(freqnames)); % Percent of time each voxel oscillates exclusively in one band, with no concurrent bands. [Sub x Vox x Band]
all_dur_secs_by_band = cell(Nsub, Nvoxin, length(freqnames)); % Stores all raw durations in seconds
all_dur_cyc_by_band  = cell(Nsub, Nvoxin, length(freqnames)); % Stores all raw durations in cycles
all_pow_by_band = cell(Nsub, Nvoxin, length(freqnames));

% Sub loop
for s=1:Nsub
    s
    load([p.data '\sub-' subs{s} '\ses-' sess{s} '\conepis\epis_v1.mat']) % to check size of time points
    Ntp = length(datasource.time);
    load([p.data '\sub-' subs{s} '\ses-' sess{s} '\episorig3cyc.mat'])


%% Transform epis from cells to matrix AND extract durations
    episodes = false(Nvoxin,Nfrex,Ntp);
    temp_dur_sec = cell(Nvoxin, length(freqnames)); 
    temp_dur_cyc = cell(Nvoxin, length(freqnames)); 
    temp_pow     = cell(Nvoxin, length(freqnames));
    
    for vx = 1:Nvoxin
        for ep = 1:length(epis3c{vx})
            fm = epis3c{vx}(ep).freq;
            fbin = dsearchn(frex',fm);
            tpts = epis3c{vx}(ep).timeps;
            episodes(vx,fbin(1),tpts) = true;
            
            % Get duration
            dur_sec = length(tpts) / srate;
            dur_cyc = dur_sec * fm;
            % Get power
            mean_pow = mean(epis3c{vx}(ep).power, 'omitmissing');

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


%% Proportion of oscillatory time 

% 1. Global
        any_tp = squeeze(sum(sum(episodes, 2), 1) > 0);
        osc_global(s) = sum(any_tp) ./ Ntp * 100;

% 2. Any frequency per voxel
        any_freq_vox = squeeze(sum(episodes, 2) > 0);
        osc_any_freq(s,:) = sum(any_freq_vox, 2) ./ Ntp * 100;

% 3. Individual frequency and voxel 
        osc_by_freq(s,:,:) = sum(episodes, 3) ./ Ntp * 100;

% 4. Frequency-Bands  
        fb_mask = false(Nvoxin, Nbands, Ntp);
        overlap_voxels = zeros(Nvoxin, Nbands, Nbands);

        for fb = 1:Nbands
            if fb == Nbands
                ct = 0; 
            else
                ct = 1;
            end
            fband = [freqbands(fb) freqbands(fb+1)];
            fbidx = dsearchn(frex', fband');
        
            fb_mask(:, fb, :) = squeeze(sum(episodes(:, fbidx(1):(fbidx(2)-ct), :), 2)) > 0;
        end


    for fb = 1:Nbands
            thisband = squeeze(fb_mask(:, fb, :));
            
            % Oscillation for each bandfreq
            osc_by_band(s, :, fb) = sum(thisband, 2) ./ Ntp * 100;
            
            % Overlaps
            for fb2 = 1:Nbands
                freq_overlap = squeeze(fb_mask(:, fb2, :));
                overlap = thisband .* freq_overlap; 
                overlap_voxels(:, fb, fb2) = sum(overlap, 2) ./ Ntp * 100;
            end
            
            band_overlap(s, :, :, :) = overlap_voxels;


            % Exclusive frequency
            idx2_fb = 1:Nbands;
            idx2_fb(fb) = []; % All bands except the current one
            other_bands = squeeze(sum(fb_mask(:, idx2_fb, :), 2)) > 0;
            
            exclusive = thisband & ~other_bands;
            
            exclusive_band(s, :, fb) = sum(exclusive, 2) ./ Ntp * 100;
        end
                
        % freq_osc_time_secs(s,:,:) = sum(episodes,3) ./ srate;
        % freq_osc_time_cycles(s,:,:) = squeeze(freq_osc_time_secs(s,:,:)) .* frex;
end

% Store results
save([p.results '\oscillatory_results_3cyc.mat'], 'osc_global', 'osc_any_freq', 'osc_by_band','osc_by_freq','band_overlap','exclusive_band', 'all_dur_secs_by_band', 'all_dur_cyc_by_band', 'all_pow_by_band')
load([p.results '\oscillatory_results_3cyc.mat'])


%% Part II. Describe and analyze

% Single value of oscillatory time in the whole brain
    mean_osc_global = mean(osc_global);
    sd_osc_global = std(osc_global, 1);

% Average of oscillatory time in each voxel at any freq
    mean_osc_any_freq = squeeze(mean(osc_any_freq, 1)); 
    sd_osc_any_freq = squeeze(std(osc_any_freq, 0, 1));

    % Figure 2b
    cfg = [];
    cfg.colmap = slanCM('BuPu');
    cfg.interp = 'linear';
    cfg.colim = [0 60];
    sBOSC_sourcefig2(mean_osc_any_freq, 'mean_osc_any_freq', cfg)
    sBOSC_nii(mean_osc_any_freq, [p.figures '\mean_osc_any_freq'])
 
     
% Frequency bands oscillatory time
    addpath('Z:\Toolbox\fieldtrip-20230118') % Path to Fieldtrip
    colmap = {'Blues', 'Greens', 'YlOrBr', 'Oranges', 'Reds'};
    Nbands = length(freqnames);
    mean_osc_band = zeros(Nvoxin, Nbands);

    for fb = 1:Nbands
            mean_osc_band(:, fb) = squeeze(mean(osc_by_band(:, :, fb), 1));
            
            cfg = [];
            cfg.colmap = slanCM(colmap{fb});
            cfg.savefig = 'no';
            cfg.interp = 'spline';
            sBOSC_sourcefig(mean_osc_band(:, fb), ['avg_' freqnames{fb} '_osc'], cfg); 
            sgtitle(['avg_' freqnames{fb} '_osc']);
            sBOSC_nii(mean_osc_band(:, fb), ['avg_' freqnames{fb} '_osc']);
    end


% Figure 2a
    osc = mean(osc_any_freq, 'all');
    nonosc = 100 - osc;
    figure, pie([nonosc osc], {'Non-oscillatory', 'Oscillatory'})
    colormap("hot")
    
    mean_per_subject = mean(osc_any_freq, 2); 
    sd_across_participants = std(mean_per_subject);
    mean_per_voxel = mean(osc_any_freq, 1);
    sd_across_voxels = std(mean_per_voxel);

    band_global_means = squeeze(mean(mean(osc_by_band, 1), 2)); 
    band_sectors = (band_global_means ./ osc) .* 100;


% Frequency band overlaps
    mean_band_overlap = squeeze(mean(band_overlap, [1, 2]));
    mean_exclusive_band = squeeze(mean(exclusive_band, [1, 2]));

% Normalize to the global oscillatory time
    overlap_graph = (mean_band_overlap / osc) * 100;
    exclusive_graph = (mean_exclusive_band / osc) * 100;

    exclusive_graph ./  diag(overlap_graph) * 100


    piecol = jet_omega_mod;
    piecols = [piecol(1,:); piecol(20,:); piecol(40,:); piecol(55,:); piecol(64,:)];
    idx_fb = 1:Nbands;
    X_coords = [0, 2, 0, -2, 0, 1.5];
    Y_coords = [0, 0, 2, 0, -2, 1.5];
    scaling_factor = 500;
    exclusive_color = [0.7 0.7 0.7]; 
    
    for fb = 1:Nbands
            figure;
            percentages = zeros(1, Nbands + 1);
            
            tmp = overlap_graph(fb, :);
            percentages(1) = tmp(fb); % Self-overlap
            
            idx2_fb = idx_fb(idx_fb ~= fb); % All other bands
            percentages(2:Nbands) = tmp(idx2_fb); % Overlap with other bands
            percentages(Nbands + 1) = exclusive_graph(fb); 
            
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
                
                if i == (Nbands + 1)
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



%% ROIs
    load(['Z:\OMEGA\OMEGA-NaturalFrequencies-main\mat_files\aal_voxel_label_10mm.mat'])
    Nroi = length(aal_label_reduc);
    
    rois = zeros(Nroi, length(freqnames), Nsub);
    roivoxs = zeros(Nvoxin, length(freqnames), Nsub);
    for roi = 1:Nroi
        for fb = 1:length(freqnames)
        inds = find(label_inside_aal_reduc==roi);
        voxs = voxel_inside_aal(inds);
        rois(roi,fb,:) = squeeze(median(osc_by_band(:,voxs,fb),2));
        roivoxs(voxs,fb,:) = repmat(squeeze(median(osc_by_band(:,voxs,fb),2))',length(voxs),1);
        end
    end
    
    % Boxplots
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
    
    sBOSC_nii(squeeze(median(roivoxs(:,3,:),3)), 'tst')

%% ROIs (osc_anyfreq)
% Table 1
    load(['Z:\OMEGA\OMEGA-NaturalFrequencies-main\mat_files\aal_voxel_label_10mm.mat'])
    Nroi = length(aal_label_reduc);

    roivoxind = zeros(1,Nvoxin); % get voxels of each roi to plot
    rois_any_freq = zeros(Nroi, Nsub);
    for roi = 1:Nroi
        inds = find(label_inside_aal_reduc==roi);
        voxs = voxel_inside_aal(inds);
        rois_any_freq(roi, :) = squeeze(median(osc_any_freq(:,voxs),2));
        roivoxind(voxs) = roi;
    end
roi_mean_any_freq = mean(rois_any_freq, 2);
roi_std_any_freq  = std(rois_any_freq, 0, 2);
roi_column = aal_label_reduc(:); 
ROI_Table_AnyFreq = table(roi_column, roi_mean_any_freq, roi_std_any_freq, ...
    'VariableNames', {'ROI_Name', 'AnyFreq_Mean_Pct', 'AnyFreq_Std_Pct'});

table_filename_anyfreq = fullfile(p.results, 'ROI_AnyFreq_Percentages.csv');
writetable(ROI_Table_AnyFreq, table_filename_anyfreq);

% Plot voxs inside ROIs
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

areas = {'vox_frontal', 'vox_parietal', 'vox_temporal', 'vox_temporal', 'vox_occipital', 'vox_medial'};

for ar = 1:length(areas)
    cfg = [];
        cfg = [];
        cfg.colmap = purplemap;
        cfg.colim = [-1 1];
        cfg.interp = 'spline';
    sBOSC_sourcefig2(eval(areas{ar}), [p.figures '/plotvoxs_' areas{ar}], cfg);
end

%% Export Raw Percentage Table for ROIs

roi_column = aal_label_reduc(:); 

ROI_Table = table(roi_column, 'VariableNames', {'ROI_Name'});

for fb = 1:Nbands
    raw_band_data = squeeze(rois(:, fb, :));
    
    band_mean = mean(raw_band_data, 2);
    band_std = std(raw_band_data, 0, 2);
    
    col_name_median = sprintf('%s_Median_Pct', freqnames{fb});
    col_name_std = sprintf('%s_Std_Pct', freqnames{fb});
    
    ROI_Table.(col_name_median) = band_mean;
    ROI_Table.(col_name_std) = band_std;
end

table_filename = fullfile(p.results, 'ROI_Raw_Percentages.csv');
writetable(ROI_Table, table_filename);

%% ROI z-scored
% Table 1
% Z-scores for each subject and band
z_osc_by_band = zeros(Nsub, Nvoxin, Nbands);

for s = 1:Nsub
    for fb = 1:Nbands
        curr_data = osc_by_band(s, :, fb);
        z_osc_by_band(s, :, fb) = (curr_data - mean(curr_data)) ./ std(curr_data);
    end
end

rois_z = zeros(Nroi, Nbands, Nsub);

for roi = 1:Nroi
    inds = find(label_inside_aal_reduc == roi);
    voxs = voxel_inside_aal(inds);
    for fb = 1:Nbands
        rois_z(roi, fb, :) = squeeze(median(z_osc_by_band(:, voxs, fb), 2));
    end
end

% Boxplots 
Ncols = 8;
Nrows = 5;
figure('Name', 'ROI z-scores by Band');
t = tiledlayout(Nrows, Ncols, 'TileSpacing', 'compact', 'Padding', 'tight');

for roi = 1:Nroi
    boxdata = squeeze(rois_z(roi, :, :))'; 
    
    nexttile(roi);
    boxplot(boxdata, 'Symbol', '', 'Widths', 0.8, 'Color', 'k', 'Notch', 'off');
    
    h = flipud(findobj(gca, 'Tag', 'Box'));
    for k = 1:length(h)
        XData = get(h(k), 'XData');
        YData = get(h(k), 'YData');    
        patch(XData, YData, piecols(k,:), 'FaceAlpha', 0.7);    
        delete(h(k));
    end
    
    hold on;
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    
    title([aal_label_reduc{roi}], 'FontSize', 8, 'Interpreter', 'none'); 
    xticklabels([]);
    set(gca, 'FontSize', 7, 'YGrid', 'on', 'GridAlpha', 0.3);
    ylim([-2.5 4]);
    box off;
end

linkaxes(t.Children, 'y');

% Another scheme
% Z-scores for each subject and band
z_osc_by_band = zeros(Nsub, Nvoxin, Nbands);
for s = 1:Nsub
    for fb = 1:Nbands
        curr_data = osc_by_band(s, :, fb);
        z_osc_by_band(s, :, fb) = (curr_data - mean(curr_data)) ./ std(curr_data);
    end
end

rois_z = zeros(Nroi, Nbands, Nsub);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc == roi);
    voxs = voxel_inside_aal(inds);
    for fb = 1:Nbands
        rois_z(roi, fb, :) = squeeze(median(z_osc_by_band(:, voxs, fb), 2));
    end
end

%% Statistics: mean of each area
region_idx   = {find(vox_frontal), find(vox_parietal), find(vox_temporal), find(vox_occipital), find(vox_medial)};
region_names = {'Frontal', 'Parietal', 'Temporal', 'Occipital', 'Medial'};
Nregions     = length(region_names);

area_data = nan(Nsub, Nbands, Nregions);

for reg = 1:Nregions
    voxs = region_idx{reg};
    for fb = 1:Nbands
        area_data(:, fb, reg) = mean(z_osc_by_band(:, voxs, fb), 2, 'omitnan');
    end
end

band_pairs  = nchoosek(1:Nbands, 2);
pair_labels = arrayfun(@(i) sprintf('%s vs %s', freqnames{band_pairs(i,1)}, freqnames{band_pairs(i,2)}), ...
              1:size(band_pairs,1), 'UniformOutput', false);

T = cell(Nregions, size(band_pairs,1));

for reg = 1:Nregions
    region_data = squeeze(area_data(:, :, reg));   
    T_rm = array2table(region_data, 'VariableNames', freqnames);

    within = table(freqnames', 'VariableNames', {'Band'});
    
    % Repeated measures
    rm  = fitrm(T_rm, 'delta-highbeta~1', 'WithinDesign', within);
    ranova_out = ranova(rm, 'WithinModel', 'Band');
    disp(region_names{reg})
    disp(ranova_out)
    
    % Post-hoc
    mc = multcompare(rm, 'Band', 'ComparisonType', 'bonferroni');
  
end

%% Boxplot
% Figure 3
figure
tl = tiledlayout(1, Nregions, 'TileSpacing', 'compact', 'Padding', 'tight');
letters = 'abcde';

for reg = 1:Nregions
    region_data = squeeze(area_data(:, :, reg)); % [Nsub x Nbands]

    % Repeated measures
    T_rm  = array2table(region_data, 'VariableNames', freqnames);
    within = table(freqnames', 'VariableNames', {'Band'});
    rm     = fitrm(T_rm, 'delta-highbeta~1', 'WithinDesign', within);
    mc     = multcompare(rm, 'Band', 'ComparisonType', 'bonferroni');

    % Multcompare: letters abcde
    letter_label = letters(1:Nbands);
    for pair = 1:size(band_pairs, 1)
        b1 = freqnames{band_pairs(pair,1)};
        b2 = freqnames{band_pairs(pair,2)};
        row = (strcmp(mc.Band_1, b1) & strcmp(mc.Band_2, b2)) | ...
              (strcmp(mc.Band_1, b2) & strcmp(mc.Band_2, b1));
        if mc.pValue(row) >= 0.05
            idx2 = band_pairs(pair,2);
            idx1 = band_pairs(pair,1);
            letter_label(idx2) = letter_label(idx1);
        end
    end

    nexttile;
    boxplot(region_data, freqnames, 'Symbol', '', 'Colors', 'k');
    h = flipud(findobj(gca, 'Tag', 'Box'));
    for k = 1:length(h)
        patch(get(h(k),'XData'), get(h(k),'YData'), piecols(k,:), 'FaceAlpha', 0.7);
        delete(h(k));
    end
    hold on;
    yline(0, '--', 'Color', [0.6 0.6 0.6]);

    y_letter = 1.5;
    for k = 1:Nbands
        text(k, y_letter, letter_label(k), ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    end

    ylim([-1.75 1.75]);
    title(region_names{reg}, 'FontSize', 11);
    if reg == 1, ylabel('Z-score'); end
    set(gca, 'FontSize', 10, 'Box', 'off', 'YGrid', 'on', 'GridAlpha', 0.3, 'TickDir', 'out');
    hold off;
end

saveas(gcf, fullfile(p.figures, 'boxplot_anova_rm_areas.png'));
close;



%% Episode duration

% Boxplots
mean_dur_cyc_subj = nan(Nsub, length(freqnames));
mean_dur_sec_subj = nan(Nsub, length(freqnames));

for s = 1:Nsub
    for fb = 1:length(freqnames)
        subj_cyc_cells = all_dur_cyc_by_band(s, :, fb);
        subj_sec_cells = all_dur_secs_by_band(s, :, fb);
        
        valid_cyc_cells = subj_cyc_cells(~cellfun(@isempty, subj_cyc_cells));
        valid_sec_cells = subj_sec_cells(~cellfun(@isempty, subj_sec_cells));
        
        if ~isempty(valid_cyc_cells)
            valid_cyc_cells = valid_cyc_cells(:);
            valid_sec_cells = valid_sec_cells(:);
            
            all_cyc_bursts = cell2mat(cellfun(@(x) x(:), valid_cyc_cells, 'UniformOutput', false));
            all_sec_bursts = cell2mat(cellfun(@(x) x(:), valid_sec_cells, 'UniformOutput', false));
            
            mean_dur_cyc_subj(s, fb) = mean(all_cyc_bursts);
            mean_dur_sec_subj(s, fb) = mean(all_sec_bursts);
        end
    end
end

grand_dur_sec_mean = nanmean(mean_dur_sec_subj, 1);
grand_dur_sec_sd = nanstd(mean_dur_sec_subj, 0, 1);
grand_dur_cyc_mean = nanmean(mean_dur_cyc_subj, 1);
grand_dur_cyc_sd = nanstd(mean_dur_cyc_subj, 0, 1);


%% rmANOVA 
% Figure 4
T_rm   = array2table(mean_dur_cyc_subj, 'VariableNames', freqnames);
within = table(freqnames', 'VariableNames', {'Band'});
rm     = fitrm(T_rm, 'delta-highbeta~1', 'WithinDesign', within);

% ANOVA
ranova_out = ranova(rm, 'WithinModel', 'Band');
disp(ranova_out)
% Post-hoc
mc = multcompare(rm, 'Band', 'ComparisonType', 'bonferroni');
t_vals = mc.Difference ./ mc.StdErr;
p_vals = mc.pValue;
df = ranova_out.DF(2); 

% Compute t
t_vals = mc.Difference ./ mc.StdErr;
p_vals = mc.pValue;

results_table = table(mc.Band_1, mc.Band_2, t_vals, p_vals, ...
    'VariableNames', {'Band1', 'Band2', 'tStat', 'pValue'});
disp(results_table)

% Boxplot
figure,
boxplot(mean_dur_cyc_subj, 'Labels', freqnames, 'Colors', 'k', 'Symbol', '.');
ylabel('Duration (Number of Cycles)');
set(gca, 'FontSize', 12, 'YGrid', 'on', 'GridAlpha', 0.3);
box off;

% Histograms
figure('Name', 'Global Distribution of Burst Durations', 'Position', [100, 100, 1200, 300], 'Color', 'w');
t = tiledlayout(1, length(freqnames), 'TileSpacing', 'compact');

for fb = 1:length(freqnames)
    band_cells = all_dur_cyc_by_band(:, :, fb);
    valid_cells = band_cells(~cellfun(@isempty, band_cells));    
    valid_cells = valid_cells(:);
    all_bursts_global = cell2mat(cellfun(@(x) x(:), valid_cells, 'UniformOutput', false));
    
    % Plot Histogram
    nexttile;
    histogram(all_bursts_global, 'BinWidth', 1, 'Normalization', 'probability', ...
        'FaceColor', piecols(fb,:), 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    
    title(freqnames{fb}, 'Interpreter', 'none', 'FontSize', 14);
    xlim([3 10]);
    ylim([0 0.6])
    xlabel('Cycles');
    
    if fb == 1
        ylabel('Probability');
    end
    set(gca, 'FontSize', 12, 'Box', 'off', 'TickDir', 'out');
end

% Seconds
figure('Name', 'Global Distribution of Burst Durations', 'Position', [100, 100, 1200, 300], 'Color', 'w');
t = tiledlayout(1, length(freqnames), 'TileSpacing', 'compact');

for fb = 1:length(freqnames)
    band_cells = all_dur_secs_by_band(:, :, fb);
    valid_cells = band_cells(~cellfun(@isempty, band_cells));
    
    valid_cells = valid_cells(:);
    all_bursts_global = cell2mat(cellfun(@(x) x(:), valid_cells, 'UniformOutput', false));
    
    nexttile;
    histogram(all_bursts_global, 'BinWidth', .2, 'Normalization', 'probability', ...
        'FaceColor', piecols(fb,:), 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    
    title(freqnames{fb}, 'Interpreter', 'none', 'FontSize', 14);
    xlim([0 2]);
    ylim([0 0.8])
    xlabel('Seconds');
    
    if fb == 1
        ylabel('Probability');
    end
    set(gca, 'FontSize', 12, 'Box', 'off', 'TickDir', 'out');
end


% Comparison
figure('Name', 'Overlaid Burst Distributions', 'Position', [100, 100, 800, 500], 'Color', 'w');
hold on;

edges = 3:1:30;
bin_centers = edges(1:end-1) + 0.5;

for fb = 1:length(freqnames)
    band_cells = all_dur_cyc_by_band(:, :, fb);
    valid_cells = band_cells(~cellfun(@isempty, band_cells));
    valid_cells = valid_cells(:);
    all_bursts_global = cell2mat(cellfun(@(x) x(:), valid_cells, 'UniformOutput', false));
    
    [N, ~] = histcounts(all_bursts_global, edges, 'Normalization', 'probability');
    
    plot(bin_centers, N, 'LineWidth', 2.5, 'Color', piecols(fb,:), 'DisplayName', freqnames{fb});
end

xlim([3 10]);
xlabel('Duration (Cycles)');
ylabel('Relative Frequency');
legend('show', 'Box', 'off');
title('Comparison of Burst Duration Distributions');
set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
hold off;

% Seconds
figure('Name', 'Overlaid Burst Distributions', 'Position', [100, 100, 800, 500], 'Color', 'w');
hold on;

edges = 0:0.2:2; % Define exact bins from 3 to 30 cycles
bin_centers = edges(1:end-1) + 0.5;

for fb = 1:length(freqnames)
    band_cells = all_dur_secs_by_band(:, :, fb);
    valid_cells = band_cells(~cellfun(@isempty, band_cells));
    valid_cells = valid_cells(:);
    all_bursts_global = cell2mat(cellfun(@(x) x(:), valid_cells, 'UniformOutput', false));
    
    % Calculate relative frequency directly using histcounts
    [N, ~] = histcounts(all_bursts_global, edges, 'Normalization', 'probability');
    
    % Plot as a thick line
    plot(bin_centers, N, 'LineWidth', 2.5, 'Color', piecols(fb,:), 'DisplayName', freqnames{fb});
end

xlim([0 2]);
xlabel('Duration (Seconds)');
ylabel('Relative Frequency');
legend('show', 'Box', 'off');
title('Comparison of Burst Duration Distributions');
set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
hold off;


%% Duration & Power 

paramlist_ep = {'dur_cyc', 'power'};

dur_subj_vox = nan(Nsub, Nvoxin, Nbands);
pow_subj_vox = nan(Nsub, Nvoxin, Nbands);

for s = 1:Nsub
    for vx = 1:Nvoxin
        for fb = 1:Nbands
            cyc_cells = all_dur_cyc_by_band{s, vx, fb};
            pow_cells  = all_pow_by_band{s, vx, fb};
            if ~isempty(cyc_cells)
                dur_subj_vox(s, vx, fb) = mean(cyc_cells, 'omitnan');
                pow_subj_vox(s, vx, fb) = mean(pow_cells, 'omitnan');
            end
        end
    end
end

cmap = slanCM('vik', 256);
cmap(100:156, :) = [];
cmap_no_white = interp1(linspace(0, 1, size(cmap,1)), cmap, linspace(0, 1, 256));

data_ep   = {dur_subj_vox, pow_subj_vox}; 

for fb = 1:Nbands
    for param = 1:length(paramlist_ep)
        tmpdata = squeeze(data_ep{param}(:, :, fb));
        zData = (tmpdata - mean(tmpdata, 2, 'omitnan')) ./ std(tmpdata, 0, 2, 'omitnan');

        % Ttest 
        [~, pvals, ~, stats] = arrayfun(@(v) ttest(zData(:,v), 0), ...
            1:Nvoxin, 'UniformOutput', false);
        pvals = cell2mat(pvals);
        tVals = cellfun(@(s) s.tstat, stats);

        % Perm
        nperm    = 1000;
        maxTval  = zeros(nperm, 1);

        for pp = 1:nperm
            permsign = (rand(Nsub, 1) > 0.5)*2 - 1;
            permdata = zData .* permsign;
            [~,~,~, permstats] = arrayfun(@(v) ttest(permdata(:,v), 0), ...
                1:Nvoxin, 'UniformOutput', false);
            permtVals  = cellfun(@(s) s.tstat, permstats);
            maxTval(pp) = max(abs(permtVals));
        end

        permTcorr    = arrayfun(@(t) (sum(maxTval >= t) + 1) / (nperm + 1), abs(tVals));
        tValspermcor = tVals;
        tValspermcor(permTcorr > 0.05) = 0;

        minval = min(tValspermcor(tValspermcor>0));
        if ~isempty(minval)
            mincolaxis(fb,param) = minval;
        end

        maxval = max(tValspermcor(tValspermcor<0));
        if ~isempty(maxval)
            maxcolaxis(fb,param) =  maxval;
        end

        % Brain map
        cfg        = [];
        cfg.colmap = cmap_no_white;
        cfg.colim  = [-10 10];
        cfg.interp = 'linear';
        cfg.savefig = 'no';
        sBOSC_sourcefig2(tValspermcor, ...
            fullfile(p.figures, [freqnames{fb} '_perm_' paramlist_ep{param}]), cfg);
        close all

    end
end

%% Correlation duration and power (Voxel-wise)
r_vox  = nan(Nvoxin, Nbands);
t_vox  = nan(Nvoxin, Nbands);
sig_vox = false(Nvoxin, Nbands);

% Inicialización de umbrales R y T
min_r_thresh = zeros(Nbands, 1);
max_r_thresh = zeros(Nbands, 1);
mincolaxis = zeros(Nbands, 1);
maxcolaxis = zeros(Nbands, 1);
cmap = flipud(slanCM('PuOr'));
cmap(100:156, :) = [];
cmap_no_white2 = interp1(linspace(0, 1, size(cmap,1)), cmap, linspace(0, 1, 256));

for fb = 1:Nbands
    dur_mat = squeeze(dur_subj_vox(:, :, fb));
    pow_mat = squeeze(pow_subj_vox(:, :, fb));
    
    dur_z = (dur_mat - mean(dur_mat, 2, 'omitnan')) ./ std(dur_mat, 0, 2, 'omitnan');
    pow_z = (pow_mat - mean(pow_mat, 2, 'omitnan')) ./ std(pow_mat, 0, 2, 'omitnan');
    
    n_valid   = sum(~isnan(dur_z) & ~isnan(pow_z), 1);
    valid_vox = n_valid >= 10;
    
    dur_c = dur_z(:, valid_vox);
    pow_c = pow_z(:, valid_vox);
    dur_c(isnan(dur_c)) = 0;
    pow_c(isnan(pow_c)) = 0;
    
    num   = sum(dur_c .* pow_c, 1);
    denom = sqrt(sum(dur_c.^2, 1) .* sum(pow_c.^2, 1));
    r_vox(valid_vox, fb) = (num ./ denom)';
    
    n = n_valid(valid_vox)';
    t_vox(valid_vox, fb) = r_vox(valid_vox, fb) .* sqrt((n-2) ./ (1 - r_vox(valid_vox, fb).^2));
    
    nperm   = 1000;
    maxTval = zeros(nperm, 1);
    for pp = 1:nperm
        permsign = (rand(Nsub, 1) > 0.5)*2 - 1;
        dur_perm = dur_c .* permsign;
        num_p    = sum(dur_perm .* pow_c, 1);
        denom_p  = sqrt(sum(dur_perm.^2, 1) .* sum(pow_c.^2, 1));
        r_perm   = (num_p ./ denom_p)';
        t_perm   = r_perm .* sqrt((n-2) ./ (1 - r_perm.^2));
        maxTval(pp) = max(abs(t_perm), [], 'omitnan');
    end
    
    t_obs       = t_vox(valid_vox, fb);
    permTcorr_r = arrayfun(@(t) (sum(maxTval >= t) + 1) / (nperm + 1), abs(t_obs));
    sig_vox(valid_vox, fb) = permTcorr_r < 0.05;
    
    r_sig_only = r_vox(:, fb);
    r_sig_only(~sig_vox(:, fb)) = 0; 
    
    min_r_val = min(r_sig_only(r_sig_only > 0));
    if ~isempty(min_r_val), min_r_thresh(fb) = min_r_val; end
    
    max_r_val = max(r_sig_only(r_sig_only < 0));
    if ~isempty(max_r_val), max_r_thresh(fb) = max_r_val; end
    
    t_vox_permcor = zeros(Nvoxin, 1);
    t_vox_permcor(valid_vox) = t_vox(valid_vox, fb);
    t_vox_permcor(~sig_vox(:, fb)) = 0;
    
    minval = min(t_vox_permcor(t_vox_permcor > 0));
    if ~isempty(minval), mincolaxis(fb) = minval; end
    
    maxval = max(t_vox_permcor(t_vox_permcor < 0));
    if ~isempty(maxval), maxcolaxis(fb) = maxval; end
    
    % --- Brain map ---
    r_masked = r_vox(:, fb);
    r_masked(~sig_vox(:, fb)) = 0;
    cfg        = [];
    cfg.colmap = cmap_no_white2;
    cfg.colim  = [-.8 .8];
    cfg.interp = 'linear';
    cfg.savefig = 'no';
    sBOSC_sourcefig2(r_masked, fullfile(p.figures, ['corr_dur_pow_' freqnames{fb}]), cfg);
    close all
end

%% Correlation duration and power (ROI)
Nroi = length(aal_label_reduc);
r_roi = nan(Nroi, Nbands);
p_roi = nan(Nroi, Nbands);

for fb = 1:Nbands
    dur_mat = squeeze(dur_subj_vox(:, :, fb));
    pow_mat = squeeze(pow_subj_vox(:, :, fb));
    
    % Z-scores
    dur_z = (dur_mat - mean(dur_mat, 2, 'omitnan')) ./ std(dur_mat, 0, 2, 'omitnan');
    pow_z = (pow_mat - mean(pow_mat, 2, 'omitnan')) ./ std(pow_mat, 0, 2, 'omitnan');
    
    for roi = 1:Nroi
        inds = find(label_inside_aal_reduc == roi);
        voxs = voxel_inside_aal(inds);
        
        x = mean(dur_z(:, voxs), 2, 'omitnan'); 
        y = mean(pow_z(:, voxs), 2, 'omitnan');
        valid = ~isnan(x) & ~isnan(y);
        
        if sum(valid) > 3
            [R_mat, P_mat] = corrcoef(x(valid), y(valid));
            r_roi(roi, fb) = R_mat(1,2);
            p_roi(roi, fb) = P_mat(1,2); 
        end
    end
end

%% Scatterplots: correlations
NCols = 5; 

for fb = 1:Nbands
    dur_mat = squeeze(dur_subj_vox(:, :, fb));
    pow_mat = squeeze(pow_subj_vox(:, :, fb));
    
    % Z-scores
    dur_z = (dur_mat - mean(dur_mat, 2, 'omitnan')) ./ std(dur_mat, 0, 2, 'omitnan');
    pow_z = (pow_mat - mean(pow_mat, 2, 'omitnan')) ./ std(pow_mat, 0, 2, 'omitnan');
    
    % Sig. rois
    sig_rois = p_roi(:, fb) < 0.05;
    r_band   = r_roi(:, fb);
    r_band(~sig_rois) = NaN;
    
    % sort
    abs_r = abs(r_band);
    [~, idx_abs] = sort(abs_r, 'descend', 'MissingPlacement', 'last');
    
    num_valid = sum(~isnan(r_band));
    roi_plot = idx_abs(1:num_valid);
       
    num_cols = NCols;
    num_rows = ceil(num_valid / num_cols);
    
    fig_width  = 300 * num_cols; 
    fig_height = 250 * num_rows;
    fig = figure('Position', [100 100 fig_width fig_height], 'Color', 'w');
    
    t = tiledlayout(num_rows, num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for ri = 1:length(roi_plot)
        roi = roi_plot(ri);
        nexttile;
        
        inds = find(label_inside_aal_reduc == roi);
        voxs = voxel_inside_aal(inds);
        x = mean(dur_z(:, voxs), 2, 'omitnan');
        y = mean(pow_z(:, voxs), 2, 'omitnan');
        valid = ~isnan(x) & ~isnan(y);
        
        col = piecols(fb,:);
        
        if p_roi(roi, fb) < 0.001
            pstr = 'p < 0.001';
        elseif p_roi(roi, fb) < 0.01
            pstr = 'p < 0.01';
        else
            pstr = 'p < 0.05'; 
        end
        
        p_fit  = polyfit(x(valid), y(valid), 1);
        x_line = linspace(min(x(valid)), max(x(valid)), 100);
        
        scatter(x(valid), y(valid), 30, col, 'filled', ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
        hold on;
        plot(x_line, polyval(p_fit, x_line), 'k-', 'LineWidth', 1.5);
        
        text(0.05, 0.90, sprintf('r = %.2f\n%s', r_roi(roi,fb), pstr), ...
            'Units', 'normalized', 'FontSize', 16, 'Color', 'k');
        
        clean_title = strrep(aal_label_reduc{roi}, '_', ' ');
        title(clean_title, 'Interpreter', 'none', 'FontSize', 20, 'FontWeight', 'normal');
        
        if ri > (num_rows - 1) * num_cols
            xlabel('Duration (Z)', 'FontSize', 10);
        end
        if mod(ri, num_cols) == 1
            ylabel('Power (Z)', 'FontSize', 10);
        end
        
        set(gca, 'FontSize', 14, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1);
        hold off;
    end
    
    save_path = fullfile(p.figures, sprintf('scatter_all_sig_rois_%s.png', freqnames{fb}));
    exportgraphics(fig, save_path, 'Resolution', 300);
    close(fig);
end