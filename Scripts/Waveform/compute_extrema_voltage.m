function [volt_peak, volt_trough] = compute_extrema_voltage(df_samples, sig)
% Compute the voltage values at the peaks and last troughs.
%
% Parameters
% ----------
% df_samples : struct
%   Struct with fields:
%       sample_peak, sample_last_trough
% sig : numeric vector
%   Time series signal
%
% Returns
% -------
% volt_peak : numeric vector
%   Voltage values at sample_peak indices
% volt_trough : numeric vector
%   Voltage values at sample_last_trough indices

    volt_peak = sig(df_samples.sample_peak);
    volt_trough = sig(df_samples.sample_last_trough);
end
