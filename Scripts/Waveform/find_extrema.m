function [peaks, troughs] = find_extrema(sig, fs, f_range, boundary, first_extrema, filter_kwargs, pass_type, pad)
% Identify peaks and troughs in a time series.
%
% MATLAB port of bycycle.cyclepoints.find_extrema (Python)
%
% Parameters
%  sig : vector
%    Time series.
%  fs : scalar
%    Sampling rate, in Hz.
%  f_range : 1x2 vector
%    Frequency range, in Hz, to narrowband filter the signal.
%  boundary : scalar, optional (default = 0)
%    Number of samples from edge of the signal to ignore.
%  first_extrema : 'peak', 'trough', or []
%    Force output to begin with a peak or trough or nothing.
%  filter_kwargs : struct, optional
%    Filter parameters. Fields: n_cycles (default 3), n_seconds (optional).
%  pass_type : string, optional (default 'bandpass')
%    Filter type.
%  pad : logical, optional (default true)
%    Whether to pad signal to prevent missed cyclepoints at edges.
%
% Returns
%  peaks : indices of peaks in sig
%  troughs : indices of troughs in sig
%

% --- INPUT PARSING & DEFAULTS ---
if nargin < 3
    error('At least sig, fs, and f_range are required.');
end
if nargin < 4 || isempty(boundary)
    boundary = 0;
end
if nargin < 5 || isempty(first_extrema)
    first_extrema = 'peak';
end
valid_first = {'peak','trough',''};
if ~ismember(first_extrema, valid_first)
    error('first_extrema must be ''peak'', ''trough'', or empty.');
end
if nargin < 6 || isempty(filter_kwargs)
    filter_kwargs = struct();
end
if ~isfield(filter_kwargs, 'n_cycles')
    filter_kwargs.n_cycles = 3;
end
if nargin < 7 || isempty(pass_type)
    pass_type = 'bandpass';
end
if nargin < 8 || isempty(pad)
    pad = true;
end

% --- Validate inputs ---
if ~isvector(sig)
    error('sig must be a vector.');
end
sig = sig(:); % ensure column vector

if ~isscalar(fs) || fs <= 0
    error('fs must be a positive scalar.');
end

if ~isvector(f_range) || numel(f_range) ~= 2 || any(f_range <= 0)
    error('f_range must be a 1x2 vector of positive values.');
end

% --- Prepare for filtering ---
sig_len = length(sig);
filt_len = 0;

if pad
    n_seconds = [];
    if isfield(filter_kwargs, 'n_seconds')
        n_seconds = filter_kwargs.n_seconds;
    end
    n_cycles = [];
    if isfield(filter_kwargs, 'n_cycles')
        n_cycles = filter_kwargs.n_cycles;
    end
    filt_len = compute_filter_length(fs, pass_type, f_range(1), f_range(2), n_cycles, n_seconds);
    pad_size = ceil(filt_len/2);
    sig = [zeros(pad_size,1); sig; zeros(pad_size,1)];
end

% --- Narrowband filter ---
b = fir1(filt_len-1, f_range / (fs/2), 'bandpass', hamming(filt_len));
sig_filt = filtfilt(b, 1, sig); % zero-phase FIR filtering

% --- Find zero-crossings ---
rise_xs = find_flank_zerox(sig_filt, 'rise');
decay_xs = find_flank_zerox(sig_filt, 'decay');

% --- Compute number of peaks and troughs ---
if rise_xs(end) > decay_xs(end)
    n_peaks = length(rise_xs) - 1;
    n_troughs = length(decay_xs);
else
    n_peaks = length(rise_xs);
    n_troughs = length(decay_xs) - 1;
end

% --- Calculate peak samples ---
peaks = zeros(n_peaks,1);
decay_xs2 = decay_xs;
for p_idx = 1:n_peaks
    last_rise = rise_xs(p_idx);
    idx = find(decay_xs2 > last_rise, 1, 'first');
    if isempty(idx)
        continue
    end
    decay_xs2 = decay_xs2(idx:end);
    next_decay = decay_xs2(1);
    [~, max_idx] = max(sig(last_rise:next_decay));
    peaks(p_idx) = max_idx - 1 + last_rise;
end

% --- Calculate trough samples ---
troughs = zeros(n_troughs,1);
rise_xs2 = rise_xs;
for t_idx = 1:n_troughs
    last_decay = decay_xs(t_idx);
    idx = find(rise_xs2 > last_decay, 1, 'first');
    if isempty(idx)
        continue
    end
    rise_xs2 = rise_xs2(idx:end);
    next_rise = rise_xs2(1);
    [~, min_idx] = min(sig(last_decay:next_rise));
    troughs(t_idx) = min_idx - 1 + last_decay;
end

% --- Remove padding ---
peaks = peaks - ceil(filt_len/2);
troughs = troughs - ceil(filt_len/2);

% --- Remove indices outside boundary ---
peaks = peaks(peaks > boundary & peaks < sig_len - boundary);
troughs = troughs(troughs > boundary & troughs < sig_len - boundary);

% --- Force first extrema if requested ---
if strcmp(first_extrema, 'peak')
    if ~isempty(troughs) && ~isempty(peaks)
        if peaks(1) > troughs(1)
            troughs(1) = [];
        end
        if peaks(end) > troughs(end)
            peaks(end) = [];
        end
    end
elseif strcmp(first_extrema, 'trough')
    if ~isempty(troughs) && ~isempty(peaks)
        if troughs(1) > peaks(1)
            peaks(1) = [];
        end
        if troughs(end) > peaks(end)
            troughs(end) = [];
        end
    end
elseif isempty(first_extrema) || strcmp(first_extrema, '')
    % do nothing
else
    error('Parameter "first_extrema" is invalid');
end

end

%% Helper function
function filt_len = compute_filter_length(fs, pass_type, f_lo, f_hi, n_cycles, n_seconds)
    if ~isempty(n_cycles) && ~isempty(n_seconds)
        error('Either `n_cycles` or `n_seconds` can be defined, but not both.');
    end

    if ~isempty(n_seconds)
        filt_len = fs * n_seconds;
    elseif ~isempty(n_cycles)
        if strcmp(pass_type, 'lowpass')
            filt_len = fs * n_cycles / f_hi;
        else
            filt_len = fs * n_cycles / f_lo;
        end
    else
        error('Either `n_cycles` or `n_seconds` needs to be defined.');
    end

    filt_len = ceil(filt_len);
    if mod(filt_len, 2) == 0
        filt_len = filt_len + 1;
    end
end
