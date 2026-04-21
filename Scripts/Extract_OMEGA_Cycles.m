%% This script applies By-Cycle algorithm adapted to MATLAB from the detected episodes with sBOSC

%% Paths

p.data = 'Z:\OMEGA\OMEGA_data';
p.results = 'Z:\Enrique\Waveform\Results';

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

load source_template_10mm_3423.mat

% Frex parameters 
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
Nfrex = length(frex);

freqbands = [frex(1) frex(9) frex(16) frex(21) frex(25) frex(30)];
freqnames = {'delta', 'theta', 'alpha', 'lowbeta', 'highbeta'};

source3423 = source;

load source_template_10mm_1925.mat

voxin = source.inverse.inside(source3423.inverse.inside);

Nvoxin = length(find(source.inverse.inside));
Ntp = 33269;

%% Apply bycyle
    ct=1;
for s = 1:Nsub
    t1 = tic;
    fprintf('\nSubject:  %d / %d\n\n', s, Nsub)
   
    load([p.data '\sub-' subs{s} '\ses-' sess{s} '\episorig3cyc.mat'])
    load([p.data '\sub-' subs{s} '\ses-' sess{s} '\datasource_3423.mat'])
    load([p.data '\sub-' subs{s} '\ses-' sess{s} '\aperiodic.mat'])

    %% Keep the 1925 voxels inside the 3423

    cfg = [];
    cfg.channel = dataclean.label(voxin);
    evalc('dataclean = ft_preprocessing(cfg, dataclean)');
    aperiodic = aperiodic(voxin,:);

    if length(dataclean.trial{1})>length(aperiodic)
        % II. Cut original signal to same length of sim_aperiodic
        cfg = [];
        cfg.begsample = 1;
        cfg.endsample = length(aperiodic);
        evalc('dataclean = ft_redefinetrial(cfg,dataclean)');
    else
        aperiodic = aperiodic(:,1:length(dataclean.trial{1}));
    end

    % III. Correction to center of the head bias
    dataclean.trial{1} = dataclean.trial{1} ./ rms(aperiodic,2);

    clear aperiodic

    fsample = dataclean.fsample;

    Ntp = size(dataclean.trial{1},2);
    time = dataclean.time{1};

    % Reconstruct episodes as matrix
    episodes = false(Nvoxin,Nfrex,Ntp); % Srate in episodes was 128, so now adapt to 256
    for vx = 1:Nvoxin
        for ep = 1:length(epis3c{vx})
            fm = epis3c{vx}(ep).freq;
            fbin = dsearchn(frex',fm);
            tpts = epis3c{vx}(ep).timeps*2;
            tpts = tpts(1):tpts(end);
            episodes(vx,fbin(1),tpts) = 1;
        end
    end

    for fb = 1:length(freqnames)

        mkdir(['Z:\Enrique\Waveform\Results\' freqnames{fb}])
        voxshape = cell(1,Nvoxin);
        cycles = cell(1,Nvoxin);

        if fb==length(freqnames) % last frequency
            ct=0;
        end

        fband = [freqbands(fb) freqbands(fb+1)];
        freq_band_idx = dsearchn(frex',fband');
        freq_band_idx = [freq_band_idx(1) freq_band_idx(end)-ct];
        fband = [frex(freq_band_idx(1)) frex(freq_band_idx(end))];

        % Voxel loop
        for vx = 1:Nvoxin
            % Organize shape features into a struct
            voxshape{vx}.time_peak = [];
            voxshape{vx}.time_trough = [];
            voxshape{vx}.volt_peak = [];
            voxshape{vx}.volt_trough = [];
            voxshape{vx}.time_decay = [];
            voxshape{vx}.time_rise = [];
            voxshape{vx}.volt_decay = [];
            voxshape{vx}.volt_rise = [];
            voxshape{vx}.sample_peak = [];
            voxshape{vx}.sample_last_zerox_decay = [];
            voxshape{vx}.sample_zerox_decay = [];
            voxshape{vx}.sample_zerox_rise = [];
            voxshape{vx}.sample_last_trough = [];
            voxshape{vx}.sample_next_trough = [];
            voxshape{vx}.volt_amp = [];
            voxshape{vx}.band_amp = [];
            voxshape{vx}.period = [];
            voxshape{vx}.time_rdsym = [];
            voxshape{vx}.time_ptsym = [];

            epis_idx = squeeze(sum(episodes(vx,freq_band_idx(1):freq_band_idx(end),:))>0);
            data = dataclean.trial{1}(vx,:);

            if ~isempty(find(epis_idx))

            % % Figure
            %  figure, plot(dataclean.time{1}, data, 'k', 'LineWidth',2)
            %  hold on
            %  plotepis = NaN(1, length(data));
            %  plotepis(epis_idx) = data(epis_idx);
            %  plot(dataclean.time{1}, plotepis,'b')
            %  xlim([23 27])
            % Bycycle

            % Filter parameters
                f_lowpass = round(mean(fband)*4); % 4 times the frequency of interest
                n_seconds_filter = .1;

            %% 1. Lowpass filter

            nyquist = fsample / 2;
            filter_len = round(n_seconds_filter * fsample);
            cutoff = f_lowpass / nyquist;

            addpath('C:\Program Files\MATLAB\R2024b\toolbox\signal\signal')
            b = fir1(filter_len, cutoff, 'low', hamming(filter_len + 1)); % lowpass FIR filter
            data_low = filtfilt(b, 1, data);

            % %Figure
            % figure, plot(dataclean.time{1}, data, 'k')
            % hold on, plot(dataclean.time{1}, data_low, 'Linewidth',3)

            %% 2. Localize peaks and troughs
            % Narrowband filter signal
            % n_seconds = .75;
            n_cycles = 3;

            [peaks, troughs] = find_extrema(data_low, fsample, fband);

            if length(peaks) == length(troughs)
            else
                continue
            end

            [rises, decays] = find_zerox(data_low, peaks, troughs);

            if length(peaks) == length(troughs) & length(peaks) == length(decays) & length(peaks) == (length(rises)+1)
            else
                continue
            end
            
            % % %Figure
            % figure, plot(time, data_low,'k')
            % hold on
            % plot(time(peaks), data_low(peaks), 'r.', 'MarkerFaceColor','r','MarkerSize',12)
            % plot(time(troughs), data_low(troughs), 'b.', 'MarkerFaceColor','b', 'MarkerSize',12)
            % plot(time(rises), data_low(rises), 'g_', 'MarkerFaceColor','g', 'MarkerSize',12)
            % plot(time(decays), data_low(decays), 'g_', 'MarkerFaceColor','g', 'MarkerSize',12)
            % 
            % Test filtering by epis: peaks and troughs

            % Peaks with epis
            peaksv = zeros(length(epis_idx),1);
            peaksv(peaks) = 1;
            peaksv2 = peaksv .* epis_idx;
            peaksepis = find(peaksv2);
            [~, locp] = ismember(peaksepis, peaks);

            % Find the cycle around peaks
            troughs_peaks = troughs(locp)';
            rises_peaks = rises(locp(1:end-1));
            decays_peaks = decays(locp);

            % Troughs with epis
            troughsv = zeros(length(epis_idx),1);
            troughsv(troughs) = 1;
            troughsv2 = troughsv .* epis_idx;
            troughsepis = find(troughsv2);
            [~, loct] = ismember(troughsepis, troughs);

            cycleidx = intersect(locp,loct);
            cycleidx = cycleidx -1;
            
            % Inflection points (2º derivative)
            inflex_idx = inflex_points2(data_low, fsample, fband);


            % For each cycle, identify the sample of each extrema and zero-crossing
            samples = [];
            samples.sample_peak = peaks(2:end)';
            samples.sample_last_zerox_decay = decays(1:end-1);
            samples.sample_zerox_decay = decays(2:end);
            samples.sample_zerox_rise = rises;
            samples.sample_last_trough = troughs(1:end-1)';
            samples.sample_next_trough = troughs(2:end)';

            num_cycles = length(samples.sample_peak);
            samples.sample_inflex_points = zeros(1, num_cycles);

            samples.sample_inflex_rise = zeros(1, num_cycles);
            samples.sample_inflex_decay = zeros(1, num_cycles);

            % Iterate over each defined cycle (Trough_i -> Peak_i -> Trough_{i+1})
            for i = 1:num_cycles
                start_search_rise = samples.sample_last_trough(i);
                end_search_rise = samples.sample_peak(i);

                inflex_rise_candidate = inflex_idx(inflex_idx > start_search_rise & inflex_idx < end_search_rise);

                if ~isempty(inflex_rise_candidate)
                    samples.sample_inflex_rise(i) = inflex_rise_candidate(end);
                else
                    samples.sample_inflex_rise(i) = samples.sample_zerox_rise(i);
                end

                start_search_decay = samples.sample_peak(i);
                end_search_decay = samples.sample_next_trough(i);

                inflex_decay_candidate = inflex_idx(inflex_idx > start_search_decay & inflex_idx < end_search_decay);

                if ~isempty(inflex_decay_candidate)
                    samples.sample_inflex_decay(i) = inflex_decay_candidate(1);
                else
                    samples.sample_inflex_decay(i) = samples.sample_last_zerox_decay(i);
                end
            end

    % Compute durations of period, peaks, and troughs
            period = samples.sample_next_trough - samples.sample_last_trough;
            time_peak = samples.sample_zerox_decay - samples.sample_zerox_rise;
            time_trough = samples.sample_zerox_rise - samples.sample_last_zerox_decay;

            % Compute extrema voltage
            [volt_peak, volt_trough] = compute_extrema_voltage(samples, data_low);

            % Compute rise-decay and peak-trough features and characteristics
            sym_features = compute_symmetry(samples, data_low, period, time_peak, time_trough);

            % Compute average oscillatory amplitude estimate during cycle
            band_amp = compute_band_amp(samples, data_low, fsample, fband, n_cycles);

            %% Take individual cycles and time
            datacycles = cell(length(cycleidx),1);
            timecyc = cell(length(cycleidx),1);
            for cp = 1:length(cycleidx)
                cyclepnts = samples.sample_last_zerox_decay(cycleidx(cp)):samples.sample_zerox_decay(cycleidx(cp));
                datacycles{cp} = data_low(cyclepnts);
                timecyc{cp} = cyclepnts;
            end

            % % Plot
            % figure
            % for ii = 1:length(datacycles)
            % plot(timecyc{ii}, datacycles{ii})
            % hold on, plot(samples.sample_peak(cycleidx(ii)),data_low(samples.sample_peak(cycleidx(ii))),'rs')
            % plot(samples.sample_last_trough(cycleidx(ii)),data_low(samples.sample_last_trough(cycleidx(ii))),'bs')
            % plot(samples.sample_last_zerox_decay(cycleidx(ii)),data_low(samples.sample_last_zerox_decay(cycleidx(ii))),'mo')
            % plot(samples.sample_zerox_decay(cycleidx(ii)),data_low(samples.sample_zerox_decay(cycleidx(ii))),'mo')
            % plot(samples.sample_zerox_rise(cycleidx(ii)),data_low(samples.sample_zerox_rise(cycleidx(ii))),'mo')
            % plot(samples.sample_zerox_rise(cycleidx(ii)),data_low(samples.sample_zerox_rise(cycleidx(ii))),'mo')
            % 
            % fprintf(['\n\nCycle: ' num2str(cycleidx(ii)) '\n']);
            % display(['Envelope: ' num2str(band_amp(cycleidx(ii)))])
            % display(['PTSYM: ' num2str(sym_features.time_ptsym(cycleidx(ii)))])
            % display(['RDSYM: ' num2str(sym_features.time_rdsym(cycleidx(ii)))])
            % 
            % pause   
            % clf
            % end


            % Store all metrics of individual cycles
            voxshape{vx}.time_peak = time_peak(cycleidx);
            voxshape{vx}.time_trough = time_trough(cycleidx);
            voxshape{vx}.volt_peak = volt_peak(cycleidx);
            voxshape{vx}.volt_trough = volt_trough(cycleidx);
            voxshape{vx}.time_decay = sym_features.time_decay(cycleidx);
            voxshape{vx}.time_rise = sym_features.time_rise(cycleidx);
            voxshape{vx}.volt_decay = sym_features.volt_decay(cycleidx);
            voxshape{vx}.volt_rise = sym_features.volt_rise(cycleidx);
            voxshape{vx}.sample_peak = samples.sample_peak(cycleidx);
            voxshape{vx}.sample_last_zerox_decay = samples.sample_last_zerox_decay(cycleidx);
            voxshape{vx}.sample_zerox_decay = samples.sample_zerox_decay(cycleidx);
            voxshape{vx}.sample_zerox_rise = samples.sample_zerox_rise(cycleidx);
            voxshape{vx}.sample_last_trough = samples.sample_last_trough(cycleidx);
            voxshape{vx}.sample_next_trough = samples.sample_next_trough(cycleidx);
            voxshape{vx}.period = period(cycleidx);
            voxshape{vx}.volt_amp = sym_features.volt_amp(cycleidx);
            voxshape{vx}.time_rdsym = sym_features.time_rdsym(cycleidx);
            voxshape{vx}.time_ptsym = sym_features.time_ptsym(cycleidx);
            voxshape{vx}.band_amp = band_amp(cycleidx);
            voxshape{vx}.sample_inflex_decay = samples.sample_inflex_decay(cycleidx);
            voxshape{vx}.sample_inflex_rise = samples.sample_inflex_rise(cycleidx);


            cycles{vx} = datacycles;
            timecycles{vx} = timecyc;

            end
        end

        % Save cycles and voxshape
        save([p.results '\' freqnames{fb} '\cycles_s' num2str(s)] , 'cycles', 'timecycles')
        save([p.results '\' freqnames{fb} '\voxshape_s' num2str(s)] , 'voxshape')

    end
            disp(['Time for sub '  num2str(s) ': ' num2str(toc(t1)/3600) ' h.'])
end

