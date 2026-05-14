function [sig_filt, kernel] = filter_signal_fir(sig, fs, pass_type, f_range, varargin)
%FILTER_SIGNAL_FIR Apply FIR filter like NeuroDSP's firwin-based filter
%
% [sig_filt, kernel] = filter_signal_fir(sig, fs, pass_type, f_range, ...)
%
% Inputs:
%   sig       - signal vector
%   fs        - sampling frequency
%   pass_type - 'bandpass', 'bandstop', 'lowpass', or 'highpass'
%   f_range   - frequency range [f_lo f_hi] or scalar cutoff
%
% Optional Name-Value:
%   'n_cycles'    - filter length in cycles (default 3)
%   'n_seconds'   - filter length in seconds (default [])
%   'remove_edges' - remove edge samples by NaN (default true)
%   'return_filter' - return filter coefficients (default false)
%
% Outputs:
%   sig_filt - filtered signal
%   kernel   - filter coefficients (if requested)

    % Parse inputs
    p = inputParser;
    addParameter(p, 'n_cycles', 3);
    addParameter(p, 'n_seconds', []);
    addParameter(p, 'remove_edges', true);
    addParameter(p, 'return_filter', false);
    parse(p, varargin{:});
    
    n_cycles = p.Results.n_cycles;
    n_seconds = p.Results.n_seconds;
    remove_edges = p.Results.remove_edges;
    return_filter = p.Results.return_filter;

    % Compute filter length like NeuroDSP
    f_lo = [];
    f_hi = [];
    if isnumeric(f_range) && length(f_range)==2
        f_lo = f_range(1);
        f_hi = f_range(2);
    elseif isnumeric(f_range) && length(f_range)==1
        if strcmp(pass_type,'lowpass')
            f_lo = [];
            f_hi = f_range;
        elseif strcmp(pass_type,'highpass')
            f_lo = f_range;
            f_hi = [];
        end
    end
    
    filt_len = compute_filter_length(fs, pass_type, f_lo, f_hi, n_cycles, n_seconds);
    
    % Design FIR filter with Kaiser window like scipy firwin defaults
    beta = 4.0; % Kaiser beta to approximate scipy's default window
    if strcmp(pass_type,'lowpass')
        Wn = f_hi/(fs/2);
        kernel = fir1(filt_len-1, Wn, 'low', kaiser(filt_len, beta), 'scale');
    elseif strcmp(pass_type,'highpass')
        Wn = f_lo/(fs/2);
        kernel = fir1(filt_len-1, Wn, 'high', kaiser(filt_len, beta), 'scale');
    elseif strcmp(pass_type,'bandpass')
        Wn = [f_lo f_hi]/(fs/2);
        kernel = fir1(filt_len-1, Wn, 'bandpass', kaiser(filt_len, beta), 'scale');
    elseif strcmp(pass_type,'bandstop')
        Wn = [f_lo f_hi]/(fs/2);
        kernel = fir1(filt_len-1, Wn, 'stop', kaiser(filt_len, beta), 'scale');
    else
        error('Unknown pass_type');
    end

    % Zero-phase filter signal
    sig_filt = filtfilt(kernel,1,sig);

    % Remove edges by NaN if requested
    if remove_edges
        edge_len = floor(filt_len/2);
        sig_filt(1:edge_len) = NaN;
        sig_filt(end-edge_len+1:end) = NaN;
    end
    
    if ~return_filter
        clear kernel
    end
end


function filt_len = compute_filter_length(fs, pass_type, f_lo, f_hi, n_cycles, n_seconds)
    % Compute filter length same as NeuroDSP
    if ~isempty(n_seconds)
        filt_len = fs * n_seconds;
    elseif ~isempty(n_cycles)
        if strcmp(pass_type,'lowpass')
            filt_len = fs * n_cycles / f_hi;
        else
            filt_len = fs * n_cycles / f_lo;
        end
    else
        error('Either n_cycles or n_seconds must be provided.');
    end
    filt_len = ceil(filt_len);
    if mod(filt_len,2)==0
        filt_len = filt_len + 1;
    end
end
