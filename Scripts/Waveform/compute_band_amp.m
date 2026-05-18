function band_amp = compute_band_amp(df_samples, sig, fs, f_range, n_cycles)
% COMPUTE_BAND_AMP  Compute mean band amplitude within each oscillatory cycle.
%
%   BAND_AMP = COMPUTE_BAND_AMP(DF_SAMPLES, SIG, FS, F_RANGE) returns the
%   mean instantaneous amplitude of the narrowband signal within each cycle
%   window, defined as the interval between consecutive troughs.
%
%   Inputs:
%       df_samples - Struct with fields:
%                       sample_peak         : peak sample indices
%                       sample_last_trough  : preceding trough per cycle
%                       sample_next_trough  : following trough per cycle
%       sig        - Time series (numeric vector)
%       fs         - Sampling rate (Hz)
%       f_range    - Frequency range [f_lo f_hi] (Hz) for bandpass filter
%       n_cycles   - Filter length in cycles at f_lo (default: 3)
%
%   Output:
%       band_amp   - Mean analytic amplitude per cycle [n_cycles x 1]

if nargin < 5 || isempty(n_cycles)
    n_cycles = 3;
end

% Compute instantaneous amplitude across the full signal.
amp = amp_by_time(sig, fs, f_range, false, n_cycles);

% Construct a trough index vector
troughs  = [df_samples.sample_last_trough(1); df_samples.sample_next_trough(:)];
n_cycles_out = length(df_samples.sample_peak);
band_amp = zeros(n_cycles_out, 1);

for c = 1:n_cycles_out
    band_amp(c) = mean(amp(troughs(c) : troughs(c + 1) -1));
end

end


%%

function amp = amp_by_time(sig, fs, f_range, remove_edges, n_cycles)
% AMP_BY_TIME  Compute instantaneous amplitude of a narrowband signal.

if nargin < 4 || isempty(remove_edges), remove_edges = true; end
if nargin < 5 || isempty(n_cycles), n_cycles = 3;    end

%  Bandpass filter
if ~isempty(f_range)
    % Filter length
    filt_len = ceil(n_cycles * fs / f_range(1));
    if mod(filt_len, 2) == 0
        filt_len = filt_len + 1;
    end
    b = fir1(filt_len - 1, f_range / (fs / 2), 'bandpass', hamming(filt_len));
    sig_filt = filtfilt(b, 1, sig);
else
    sig_filt = sig;
    filt_len = 1;
end

%  Instantaneous amplitude via Hilbert transform
amp = abs(hilbert(sig_filt));

if remove_edges && filt_len > 1
    half_len          = floor(filt_len / 2);
    amp(1:half_len)   = NaN;
    amp(end - half_len + 1 : end) = NaN;
end

end