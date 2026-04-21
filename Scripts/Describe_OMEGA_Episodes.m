%% This script loads previously identified episodes with sBOSC with a duration >= 3 cycles and extracts descriptives stats and figures.

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

Nvoxin = length(find(source.inverse.inside));
Ntp = 33269;

%% Part I: Extract data from episodes
osc_global = zeros(Nsub,1); % Percent of total time where at least one voxel is oscillating at any frequency [Sub x 1]
osc_any_freq = zeros(Nsub,Nvoxin); % Percent of time each voxel is oscillating (collapsed across all frequencies). [Sub x Vox]
osc_by_freq = zeros(Nsub,Nvoxin, length(frex)-2); % Percent of time each voxel is oscillating at each frequency. [Sub x Vox x Freq]
osc_by_band = zeros(Nsub,Nvoxin, 5); % Percent of time each voxel is oscillating within frequency bands. [Sub x Vox x Band]
band_overlap =  zeros(Nsub,Nvoxin,length(freqnames),length(freqnames)); % Percent of temporal overlap in oscillations between pairs of bands per voxel. [Sub x Vox x Band1 x Band2]
exclusive_band = zeros(Nsub,Nvoxin, length(freqnames)); % Percent of time each voxel oscillates exclusively in one band, with no concurrent bands. [Sub x Vox x Band]
 
% Sub loop
for s=1:Nsub
    s
    load([p.data '\sub-' subs{s} '\ses-' sess{s} '\episorig3cyc.mat'])

%% Transform epis from cells to matrix
    episodes = false(Nvoxin,Nfrex,Ntp); 
    for vx = 1:Nvoxin
        for ep = 1:length(epis3c{vx})
            fm = epis3c{vx}(ep).freq;
            fbin = dsearchn(frex',fm);
            tpts = epis3c{vx}(ep).timeps;
            episodes(vx,fbin(1),tpts) = true;
        end
    end

    % Cut at 33 Hz
    episodes = episodes(:,1:30,:);

%% Proportion of oscillatory time 
Ntp = size(episodes, 3);
Nbands = length(freqnames);

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
save([p.results '\oscillatory_results_3cycb', 'osc_global', 'osc_any_freq', 'osc_by_band','osc_by_freq','band_overlap','exclusive_band'])

load([p.results '\oscillatory_results_3cycb.mat'])


%% Part II. Describe and analyze

% Single value of oscillatory time in the whole brain
    mean_osc_global = mean(osc_global);
    sd_osc_global = std(osc_global, 1);

% Average of oscillatory time in each voxel at any freq
    mean_osc_any_freq = squeeze(median(osc_any_freq, 1)); 
    sd_osc_any_freq = squeeze(std(osc_any_freq, 0, 1));

    % Get .nii figures
    sBOSC_nii(mean_osc_any_freq, 'mean_osc_any_freq')
    sBOSC_nii(sd_osc_any_freq, 'sd_osc_any_freq')
     
    % Plot a figure
    cfg = [];
    cfg.colmap = slanCM('speed');
    cfg.interp = 'linear';
    cfg.colim = [0 100];
    sBOSC_sourcefig(mean_osc_any_freq, 'mean_osc_any_freq', cfg) % add if filename save

% Frequency bands oscillatory time
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


% Percentages graph
    osc = mean(mean_osc_any_freq);
    nonosc = 100 - osc;
    figure, pie([nonosc osc], {'Non-oscillatory', 'Oscillatory'})
    colormap("hot")

    band_global_means = squeeze(mean(mean(osc_by_band, 1), 2)); 
    band_sectors = (band_global_means ./ osc) .* 100;


% Frequency band overlaps
    mean_band_overlap = squeeze(mean(band_overlap, [1, 2]));
    mean_exclusive_band = squeeze(mean(exclusive_band, [1, 2]));

% Normalize to the global oscillatory time
    overlap_graph = (mean_band_overlap / osc) * 100;
    exclusive_graph = (mean_exclusive_band / osc) * 100;

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
