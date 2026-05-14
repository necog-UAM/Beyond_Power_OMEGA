function band_amp = compute_band_amp(df_samples, sig, fs, f_range, n_cycles)
    % compute_band_amp Compute the average amplitude of each oscillation.
    %
    % Parameters
    % ----------
    % sig : 1d array
    %     Time series.
    % fs : float
    %     Sampling rate, in Hz.
    % f_range : 2-element vector
    %     Frequency range for narrowband signal of interest (Hz).
    % n_cycles : int, optional, default: 3
    %     Length of filter, in number of cycles, at the lower cutoff frequency.
    %
    % Returns
    % -------
    % band_amp : 1d array
    %     Average analytic amplitude of the oscillation.
    %
    % Example:
    % fs = 500;
    % sig = sim_bursty_oscillation(10, fs, 10); % Assuming similar function exists
    % df_samples = compute_cyclepoints(sig, fs, [8 12]);
    % band_amp = compute_band_amp(df_samples, sig, fs, [8 12]);

    if nargin < 5
        n_cycles = 3;
    end

    amp = amp_by_time(sig, fs, f_range, false, n_cycles);


    % Construct troughs vector
    troughs = [df_samples.sample_last_trough(1); df_samples.sample_next_trough];

    band_amp = zeros(length(df_samples.sample_peak), 1);
    for sig_idx = 1:length(df_samples.sample_peak)
        band_amp(sig_idx) = mean(amp(troughs(sig_idx) : troughs(sig_idx + 1)));
    end

    % figure, plot(sig(df_samples.sample_last_trough(sig_idx):df_samples.sample_next_trough(sig_idx)))
    % hold on
    % plot(amp(df_samples.sample_last_trough(sig_idx):df_samples.sample_next_trough(sig_idx)))
    % 

end
