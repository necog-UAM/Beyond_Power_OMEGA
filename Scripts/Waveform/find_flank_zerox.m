function zero_xs = find_flank_zerox(sig, flank, midpoint)
% FIND_FLANK_ZEROX Find zero-crossings on rising or decaying flanks of a signal.
%
% Parameters
% ----------
% sig : 1D array
%     Time series to detect zero-crossings in.
% flank : 'rise' or 'decay'
%     Which flank, 'rise' or 'decay', to use to find zero-crossings.
% midpoint : scalar (optional)
%     Midpoint value to compare against. Default is 0.
%
% Returns
% -------
% zero_xs : 1D array
%     Indices of zero-crossings in the signal.

    if nargin < 3 || isempty(midpoint)
        midpoint = 0;
    end

    if ~strcmp(flank, 'rise') && ~strcmp(flank, 'decay')
        error('flank must be ''rise'' or ''decay''');
    end

    if strcmp(flank, 'rise')
        pos = sig <= midpoint;
    else
        pos = sig > midpoint;
    end

    % Find zero-crossing indices
    zero_xs = find(pos(1:end-1) & ~pos(2:end));

    % If no zero-crossings found, return middle index
    if isempty(zero_xs)
        zero_xs = floor(length(sig) / 2);
    end
end
