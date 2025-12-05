function [dtw_fii, tr_dtw] = compute_dtw_fii(x, y, gamma, FWHM, win_size)
% compute_dtw_fii - Computes FII and time-resolved DTW for a pair of signals.
%
% Syntax:
%   [dtw_fii, tr_dtw] = compute_group_dtw(x, y, gamma, Tr, band)
%
% Inputs:
%   x         - A single time course
%   y         - A single time course
%   gamma     - The gamma parameter to use in the DTW computation.
%   FWHM      - Full width half maximum for Gaussian smoothing.
%   win_size  - Window size for DTW (in samples) and anchor for FII.
%
% Outputs:
%   dtw_fii   - Functional inertia index values.
%   tr_dtw    - Time-resolved DTW distances
%
% Description:
%   This function computes the functional inertia index from the time-resolved DTW distances. 
%   The DTWFramework function is used to compute time-resolved DTW.
%
% Example:
%   % Assume x and y are two time series,
%   % gamma = 1.5; FWHM = 5; win_size = 4;
%   [full_dtw, full_dtw_norm] = compute_dtw_fii(x, y, gamma, FWHM, win_size);

    % Determine dimensions from input
    num_timepoints = max([length(x) length(y)]);

    % choose anchor index in samples (e.g. based on window size)
    t_a = win_size;
    num_cumm = num_timepoints - t_a + 1;

    % Display progress information in the command window.
    fprintf('Processing gamma = %.2f \n', gamma);
        
    % Normalize the time courses.
    x = zscore(x);
    y = zscore(y);
    
    % Compute time-resolved DTW using the custom DTW framework.
    [~, ~, tr_dtw, ~] = DTWFramework(x, y, gamma, win_size);

    dtw_cumm = zeros(num_cumm, 1);
    
    for t = t_a:num_timepoints
        dtw_cumm(t - t_a + 1) = mean(tr_dtw(1:t));
    end
    
    smoothed_cumm = gauss_smooth(dtw_cumm, FWHM);
    dtw_fii = diff(smoothed_cumm);

end
