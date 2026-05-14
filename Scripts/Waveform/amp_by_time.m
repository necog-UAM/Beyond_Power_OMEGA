function amp = amp_by_time(sig, fs, f_range, remove_edges, n_cycles)
    % Compute instantaneous amplitude of a time series
    % Inputs:
    %   sig          - input signal (vector)
    %   fs           - sampling frequency
    %   f_range      - [low, high] frequency range for bandpass
    %   remove_edges - boolean, if true set edge samples to NaN
    %   n_cycles     - filter length in number of cycles at low freq
    
    if nargin < 4
        remove_edges = true;
    end
    if nargin < 5
        n_cycles = 3; 
    end
    
    if ~isempty(f_range)
        % Design FIR bandpass filter
        filt_len = round(n_cycles * fs / f_range(1)); % filter length = n_cycles at low freq
        if mod(filt_len,2) == 0
            filt_len = filt_len + 1; % make filter length odd for fir1
        end
        
        b = fir1(filt_len-1, f_range/(fs/2), hamming(filt_len)); % bandpass FIR filter
        % Zero-phase filtering
        sig_filt = filtfilt(b, 1, sig);
    else
        sig_filt = sig; % no filtering
        filt_len = 1;
    end
    
    % Analytic signal via Hilbert transform
    analytic_sig = hilbert(sig_filt);
    
    % Instantaneous amplitude = abs of analytic signal
    amp = abs(analytic_sig);
    
    % Remove edges if requested
    if remove_edges && filt_len > 1
        half_len = floor(filt_len/2);
        amp(1:half_len) = NaN;
        amp(end-half_len+1:end) = NaN;
    end
end
