function sym_features = compute_symmetry(df_samples, sig, period, time_peak, time_trough)
% COMPUTE_SYMMETRY  Compute rise-decay and peak-trough symmetry features.
%
%   SYM_FEATURES = COMPUTE_SYMMETRY(DF_SAMPLES, SIG) computes waveform
%   symmetry metrics: period, time_peak, and time_trough.
%
%   Inputs:
%       df_samples  - Struct with fields:
%                       sample_peak, sample_last_trough, sample_next_trough,
%                       sample_zerox_rise, sample_zerox_decay,
%                       sample_last_zerox_decay
%       sig         - Filtered time series (numeric vector)
%       period      - Cycle period in samples [optional]
%       time_peak   - Peak half-cycle duration in samples [optional]
%       time_trough - Trough half-cycle duration in samples [optional]
%
%   Output:
%       sym_features - Struct with fields:
%                       time_rise, time_decay  : half-cycle durations (samples)
%                       volt_rise, volt_decay  : voltage excursions from peak
%                       volt_amp               : mean peak-to-trough amplitude
%                       time_rdsym             : rise-decay symmetry (rise / period)
%                       time_ptsym             : peak-trough symmetry
%                                                (time_peak / (time_peak + time_trough))


if nargin < 5 || isempty(period) || isempty(time_peak) || isempty(time_trough)
    period = df_samples.sample_next_trough - df_samples.sample_last_trough;
    time_peak = df_samples.sample_zerox_decay - df_samples.sample_zerox_rise;
    time_trough = df_samples.sample_zerox_rise - df_samples.sample_last_zerox_decay;
end


%  Rise and decay durations (samples)
time_rise = df_samples.sample_peak - df_samples.sample_last_trough;
time_decay = df_samples.sample_next_trough - df_samples.sample_peak;


volt_rise = sig(df_samples.sample_peak) - sig(df_samples.sample_last_trough);
volt_decay = sig(df_samples.sample_peak) - sig(df_samples.sample_next_trough);
volt_amp = (volt_rise + volt_decay) / 2;


%  Symmetry metrics
time_rdsym = time_rise  ./ period;
time_ptsym = time_peak  ./ (time_peak + time_trough);

sym_features = struct( ...
    'time_rise',  time_rise,  ...
    'time_decay', time_decay, ...
    'volt_rise',  volt_rise,  ...
    'volt_decay', volt_decay, ...
    'volt_amp',   volt_amp,   ...
    'time_rdsym', time_rdsym, ...
    'time_ptsym', time_ptsym);

end