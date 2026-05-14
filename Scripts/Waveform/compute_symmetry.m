function sym_features = compute_symmetry(df_samples, sig, period, time_peak, time_trough)
% Compute rise-decay and peak-trough symmetry features from cycle points.
%
% Parameters
% ----------
% df_samples : struct
%   Struct with fields:
%       sample_peak, sample_last_trough, sample_next_trough
% sig : 1D numeric array
%   Time series signal.
% period : optional, numeric array
%   Periods of each cycle.
% time_peak : optional, numeric array
%   Time from rise to decay.
% time_trough : optional, numeric array
%   Time from previous trough to rise.
%
% Returns
% -------
% sym_features : struct
%   Contains symmetry features:
%       time_decay, time_rise, volt_decay, volt_rise, volt_amp,
%       time_rdsym, time_ptsym

    % Compute rise and decay durations
    time_decay = df_samples.sample_next_trough - df_samples.sample_peak;
    time_rise = df_samples.sample_peak - df_samples.sample_last_trough;

    % Compute voltage differences
    volt_decay = sig(df_samples.sample_peak) - sig(df_samples.sample_next_trough);
    volt_rise = sig(df_samples.sample_peak) - sig(df_samples.sample_last_trough);
    volt_amp = (volt_decay + volt_rise) / 2;

    % Compute durations if not provided
    if nargin < 3 || isempty(period) || isempty(time_peak) || isempty(time_trough)   
    period = df_samples.sample_next_trough - df_samples.sample_last_trough;
    time_peak = df_samples.sample_zerox_decay - df_samples.sample_zerox_rise;
    time_trough = df_samples.sample_zerox_rise - df_samples.sample_last_zerox_decay;
    end

    % Compute symmetry metrics
    time_rdsym = time_rise ./ period;
    time_ptsym = time_peak ./ (time_peak + time_trough);

    % Store features in struct
    sym_features.time_decay = time_decay;
    sym_features.time_rise = time_rise;
    sym_features.volt_decay = volt_decay;
    sym_features.volt_rise = volt_rise;
    sym_features.volt_amp = volt_amp;
    sym_features.time_rdsym = time_rdsym;
    sym_features.time_ptsym = time_ptsym;
end
