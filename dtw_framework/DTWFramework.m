function [dtw_d, dtw_d_norm, tr_dtw, d_dtw] = DTWFramework(x, y, gamma, win_size, num_timepoints)
% DTWFramework - Computes DTW distance and extracts time-resolved amplitude mismatch
%                between two signals after DTW warping.
%
% Syntax:
%   [dtw_d, dtw_d_norm, tr_dtw, d_dtw] = DTWFramework(x, y, gamma, win_size, num_timepoints)
%
% Inputs:
%   x, y          - Input time course vectors.
%   gamma         - Exponent used in DTW computations.
%   win_size      - Window size used for DTW computation.
%   num_timepoints- (Optional) Desired number of time points for interpolation. If not
%                   provided or empty, defaults to the maximum of length(x) and length(y).
%
% Outputs:
%   dtw_d       - DTW distance between x and y.
%   dtw_d_norm  - DTW distance normalized by the number of alignment steps.
%   tr_dtw      - Time-resolved DTW (absolute interpolated difference of amplitude mismatch).
%   d_dtw       - Directional DTW (difference tr_DTW reference x -
%   reference y)
%
% Description:
%   This function calculates the DTW distance between two signals and extracts the
%   time-resolved amplitude mismatch after DTW warping. If the desired number of
%   time points for interpolation is not specified, it defaults to the longer signal.
%
% Example:
%   [dtw_d, dtw_d_norm, tr_dtw, d_dtw] = DTWFramework(x, y, 1, 20);

    % Set defaults if not provided or empty.
    if nargin < 4 || isempty(win_size)
        win_size = max(length(x), length(y));
    end

    if nargin < 5 || isempty(num_timepoints)
        num_timepoints = max(length(x), length(y));
    end

    
    % Compute the DTW distance and alignment indices using a custom DTW function.
    [~, dtw_d, ix, iy] = dtw_custom(x, y, win_size, gamma);
    dtw_d_norm = dtw_d / length(ix);
    
    % Extract time-resolved amplitude mismatch and reference vectors.
    [tr_dtw, d_dtw_x, d_dtw_y] = extract_tr_DTW(x, y, ix, iy, num_timepoints, gamma);
    
    % Compute the difference between the reference vectors.
    d_dtw = d_dtw_x - d_dtw_y;
end
