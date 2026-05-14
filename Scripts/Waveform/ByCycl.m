% Cycle by cycle approach
load C:\Users\Usuario\Desktop\sig.mat
% sig = dataclean.trial{1}(vx,:);
% fsample = 128;
fsample=1250;
times = 0:1/fsample:(length(sig)/fsample - 1/fsample);
% Filter parameters
f_theta = [8, 12]
f_lowpass = 40
n_seconds_filter = .1

%% 1. Lowpass filter

nyquist = fsample / 2;
filter_len = round(n_seconds_filter * fsample);  % filter length in samples
cutoff = f_lowpass / nyquist;

b = fir1(filter_len, cutoff, 'low', hamming(filter_len + 1)); % lowpass FIR filter
sig_low = filtfilt(b, 1, sig); % zero-phase filtering

%% 2. Localize peaks and troughs
% Narrowband filter signal
n_seconds_theta = .75
n_cycles = 3;

[peaks, troughs] = find_extrema(sig_low, fsample, f_theta);
[rises, decays] = find_zerox(sig_low, peaks, troughs);

figure, plot(times, sig_low)
hold on
plot(times(peaks), sig_low(peaks), 'ro', 'MarkerFaceColor','r')
plot(times(troughs), sig_low(troughs), 'bo', 'MarkerFaceColor','b')
plot(times(rises), sig_low(rises), 'yo', 'MarkerFaceColor','y')
plot(times(decays), sig_low(decays), 'go', 'MarkerFaceColor','g')

%% For each cycle, identify the sample of each extrema and zero-crossing
samples = [];
samples.sample_peak = peaks(2:end);                  
samples.sample_last_zerox_decay = decays(1:end-1);  
samples.sample_zerox_decay = decays(2:end);         
samples.sample_zerox_rise = rises;                    
samples.sample_last_trough = troughs(1:end-1);       
samples.sample_next_trough = troughs(2:end);     

% Compute durations of period, peaks, and troughs
[period, time_peak, time_trough] = compute_durations(samples);

% Compute extrema voltage
[volt_peak, volt_trough] = compute_extrema_voltage(samples, sig);

% Compute rise-decay and peak-trough features and characteristics
sym_features = compute_symmetry(samples, sig, period, time_peak, time_trough);

% Compute average oscillatory amplitude estimate during cycle
band_amp = compute_band_amp(samples, sig, fsample, f_theta, n_cycles); 

% Organize shape features into a struct 
shape.period = period;
shape.time_peak = time_peak;
shape.time_trough = time_trough;
shape.volt_peak = volt_peak;
shape.volt_trough = volt_trough;
shape.time_decay = sym_features.time_decay;
shape.time_rise = sym_features.time_rise;
shape.volt_decay = sym_features.volt_decay;
shape.volt_rise = sym_features.volt_rise;
shape.volt_amp = sym_features.volt_amp;
shape.time_rdsym = sym_features.time_rdsym;
shape.time_ptsym = sym_features.time_ptsym;
shape.band_amp = band_amp;
shape.sample_peak = samples.sample_peak;               
shape.sample_last_zerox_decay = samples.sample_last_zerox_decay;
shape.sample_zerox_decay = samples.sample_zerox_decay;
shape.sample_zerox_rise = samples.sample_zerox_rise;
shape.sample_last_trough = samples.sample_last_trough;
shape.sample_next_trough = samples.sample_next_trough;

%% Filter by epis
