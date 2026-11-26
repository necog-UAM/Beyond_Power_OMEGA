function [peaks, troughs] = find_extrema(sig, fs, f_range, boundary, first_extrema, filter_opts, pass_type, pad)
% FIND_EXTREMA  Identify peaks and troughs in a time series.
%
%   [PEAKS, TROUGHS] = FIND_EXTREMA(SIG, FS, F_RANGE) returns sample
%   indices of peaks and troughs in SIG after narrowband filtering within
%   F_RANGE.
%
%   Inputs:
%       sig           - Time series (vector)
%       fs            - Sampling rate (Hz)
%       f_range       - Frequency range [f_lo f_hi] (Hz)
%       boundary      - Samples to ignore at signal edges (default: 0)
%       first_extrema - Force output to begin with 'peak', 'trough', or ''
%                       (default: 'peak')
%       filter_opts - Struct with optional fields:
%                           n_cycles  (default: 3)
%                           n_seconds (overrides n_cycles if set)
%       pass_type     - Filter type string (default: 'bandpass')
%       pad           - Logical; pad signal edges to reduce edge artefacts
%                       (default: true)
%
%   Outputs:
%       peaks   - Sample indices of peaks in sig
%       troughs - Sample indices of troughs in sig


%  Input validation and defaults
if nargin < 3
    error('find_extrema:insufficientInputs', 'At least sig, fs, and f_range are required.');
end
if nargin < 4 || isempty(boundary),      boundary      = 0;          end
if nargin < 5 || isempty(first_extrema), first_extrema = 'peak';     end
if nargin < 6 || isempty(filter_opts), filter_opts = struct();    end
if nargin < 7 || isempty(pass_type),     pass_type     = 'bandpass'; end
if nargin < 8 || isempty(pad),           pad           = true;        end

if ~ismember(first_extrema, {'peak', 'trough', ''})
    error('find_extrema:invalidInput', 'first_extrema must be ''peak'', ''trough'', or ''''.');
end
if ~isvector(sig)
    error('find_extrema:invalidInput', 'sig must be a vector.');
end
if ~isscalar(fs) || fs <= 0
    error('find_extrema:invalidInput', 'fs must be a positive scalar.');
end
if ~isvector(f_range) || numel(f_range) ~= 2 || any(f_range <= 0)
    error('find_extrema:invalidInput', 'f_range must be a 2-element vector of positive values.');
end
if ~isfield(filter_opts, 'n_cycles'),  filter_opts.n_cycles  = 3;   end
if ~isfield(filter_opts, 'n_seconds'), filter_opts.n_seconds = [];  end

sig     = sig(:);   
sig_len = length(sig);

%  Compute filter length 
filt_len = compute_filter_length(fs, pass_type, f_range(1), f_range(2), filter_opts.n_cycles, filter_opts.n_seconds);

%  Pad signal
if pad
    pad_size = ceil(filt_len / 2);
    sig      = [zeros(pad_size, 1); sig; zeros(pad_size, 1)];
end

%  Narrowband filter
b        = fir1(filt_len - 1, f_range / (fs / 2), 'bandpass', hamming(filt_len));
sig_filt = filtfilt(b, 1, sig);


%  Find zero-crossings on rising and decaying flanks
rise_xs  = find_flank_zerox(sig_filt, 'rise');
decay_xs = find_flank_zerox(sig_filt, 'decay');

%  Determine number of peaks and troughs
if rise_xs(end) > decay_xs(end)
    n_peaks   = length(rise_xs) - 1;
    n_troughs = length(decay_xs);
else
    n_peaks   = length(rise_xs);
    n_troughs = length(decay_xs) - 1;
end

%  Locate peak samples
peaks = zeros(n_peaks, 1);
for p_idx = 1:n_peaks
    last_rise      = rise_xs(p_idx);
    next_decay_idx = find(decay_xs > last_rise, 1, 'first');
    if isempty(next_decay_idx), continue; end
    next_decay     = decay_xs(next_decay_idx);
    [~, max_idx]   = max(sig(last_rise:next_decay-1));
    peaks(p_idx)   = last_rise + max_idx - 1;
end

%  Locate trough samples
troughs = zeros(n_troughs, 1);
for t_idx = 1:n_troughs
    last_decay    = decay_xs(t_idx);
    next_rise_idx = find(rise_xs > last_decay, 1, 'first');
    if isempty(next_rise_idx), continue; end
    next_rise     = rise_xs(next_rise_idx);
    [~, min_idx]  = min(sig(last_decay:next_rise-1));
    troughs(t_idx) = last_decay + min_idx - 1;
end

%  Remove padding offset and boundary samples
if pad
    peaks = peaks - ceil(filt_len / 2);
    troughs = troughs - ceil(filt_len / 2);
end
peaks   = peaks(peaks   > boundary & peaks   <= sig_len - boundary);
troughs = troughs(troughs > boundary & troughs <= sig_len - boundary);


%  Force first extrema type
if strcmp(first_extrema, 'peak')
    if ~isempty(peaks) && ~isempty(troughs)
        if peaks(1)   > troughs(1),   troughs(1)   = []; end
        if peaks(end) > troughs(end), peaks(end)   = []; end
    end
elseif strcmp(first_extrema, 'trough')
    if ~isempty(peaks) && ~isempty(troughs)
        if troughs(1)   > peaks(1),   peaks(1)     = []; end
        if troughs(end) > peaks(end), troughs(end) = []; end
    end
end
end



%% Local functions
function filt_len = compute_filter_length(fs, pass_type, f_lo, f_hi, n_cycles, n_seconds)
if ~isempty(n_cycles) && ~isempty(n_seconds)
    error('compute_filter_length:ambiguousInput', 'Specify either n_cycles or n_seconds.');
end
if ~isempty(n_seconds)
    filt_len = round(fs * n_seconds);
elseif ~isempty(n_cycles)
    if strcmp(pass_type, 'lowpass')
        filt_len = ceil(fs * n_cycles / f_hi);
    else
        filt_len = ceil(fs * n_cycles / f_lo);
    end
else
    error('compute_filter_length:missingInput', 'Either n_cycles or n_seconds must be specified.');
end
% Force odd length
if mod(filt_len, 2) == 0
    filt_len = filt_len + 1;
end
end



%%

function zero_xs = find_flank_zerox(sig, flank, midpoint)

if nargin < 3 || isempty(midpoint)
    midpoint = 0;
end
if ~ismember(flank, {'rise', 'decay'})
    error('find_flank_zerox:invalidInput', ...
        'flank must be ''rise'' or ''decay''.');
end
if strcmp(flank, 'rise')
    pos = sig <= midpoint;   
else
    pos = sig >  midpoint;  
end

zero_xs = find(pos(1:end-1) & ~pos(2:end));

if isempty(zero_xs)
    zero_xs = floor(length(sig) / 2) + 1;
end

end