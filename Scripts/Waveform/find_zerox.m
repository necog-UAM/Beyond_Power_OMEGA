function [rises, decays] = find_zerox(sig, peaks, troughs)
    n_rises = length(peaks);
    n_decays = length(troughs);
    idx_bias = 0;

    if peaks(1) < troughs(1)
        n_rises = n_rises - 1;
    else
        n_decays = n_decays - 1;
        idx_bias = idx_bias + 1;
    end

    rises = find_flank_midpoints(sig, 'rise', n_rises, troughs, peaks, idx_bias);
    decays = find_flank_midpoints(sig, 'decay', n_decays, peaks, troughs, idx_bias);
end

%% Local function
function flanks = find_flank_midpoints(sig, flank, n_flanks, extrema_start, extrema_end, idx_bias)
    if strcmp(flank, 'rise')
        idx_bias = -idx_bias + 1;
        comp = @(x, y) x > y;
    else
        comp = @(x, y) x < y;
    end

    flanks = zeros(1, n_flanks);
    for idx = 1:n_flanks
        i_start = extrema_start(idx);
        i_end = extrema_end(idx + idx_bias);

        sig_temp = sig(i_start:i_end);

        if sum(abs(sig_temp)) == 0 || comp(sig_temp(1), sig_temp(end))
            flanks(idx) = i_start + floor(length(sig_temp)/2);
        else
            midpoint = (sig_temp(1) + sig_temp(end)) / 2;
            zc = find_flank_zerox(sig_temp, flank, midpoint);
            flanks(idx) = i_start + floor(median(zc)) -1;
        end
    end
end

function zero_xs = find_flank_zerox(sig, flank, midpoint)
    if nargin < 3
        midpoint = 0;
    end

    if strcmp(flank, 'rise')
        pos = sig <= midpoint;
    elseif strcmp(flank, 'decay')
        pos = sig > midpoint;
    else
        error('Flank must be ''rise'' or ''decay''.');
    end

    zero_xs = find(pos(1:end-1) & ~pos(2:end));

    if isempty(zero_xs)
        zero_xs = floor(length(sig)/2) +1;
    end
end
