function [full_dtw_fii, full_tr_dtw] = compute_subject_dtw_fii(subTcs, gamma, FWHM, win_size)
% compute_subject_dtw_fii - Computes FII and time-resolved DTW metrics for all subjects and unique component pairs.
%
% Syntax:
%   [full_dtw_fii, full_tr_dtw] = compute_group_dtw(subTcs, gamma, Tr, band)
%
% Inputs:
%   subTcs    - A 3D matrix of time courses with dimensions 
%               [num_subjects x num_timepoints x num_components].
%   gamma     - The gamma parameter to use in the DTW computation.
%   FWHM      - Full width half maximum for Gaussian smoothing.
%   win_size  - Window size for DTW (in samples) and anchor for FII.
%
% Outputs:
%   full_dtw_fii    - A 3D array of functional inertia index values with
%   dimensions [num_subjects x num_timepoints-win_size x num_features].
%                     [num_gammas x num_subjects x num_features], where num_features = nchoosek(num_components, 2).
%   full_tr_dtw     - A 3D array of time-resolved DTW distances with
%   dimensions [num_subjects x num_timepoints x num_features]..
%
% Description:
%   This function computes the functional inertia index from the time-resolved DTW distances 
%   between all unique pairs of components for each subject. For each subject, it computes a 
%   lower-triangular matrix of symmetric DTW distances between the components, vectorizes the 
%   unique values using icatb_mat2vec, and stores the result. The DTWFramework function is used 
%   to compute time-resolved DTW.
%
% Example:
%   % Assume subTcs is a 3D matrix of dimensions [100 subjects x 150 timepoints x 10 components],
%   % gamma = 1.5; FWHM = 5; win_size = 4;
%   [full_dtw, full_dtw_norm] = compute_subject_dtw_fii(subTcs, gamma, FWHM, win_size);

    % Determine dimensions from input
    num_subjects   = size(subTcs, 1);
    num_timepoints = size(subTcs, 2);
    num_components = size(subTcs, 3);
    num_features   = nchoosek(num_components, 2);

    % Preallocate output matrices
    full_tr_dtw = zeros(num_subjects, num_timepoints, num_features);
    full_dtw_fii = zeros(num_subjects, num_timepoints-win_size, num_features);
    
    % choose anchor index in samples (e.g. based on window size)
    t_a = win_size;
    num_cumm = num_timepoints - t_a + 1;

    % Loop over each subject
    for sub = 1:num_subjects
        % Display progress information in the command window.
        fprintf('Processing gamma = %.2f, subject = %d/%d\n', ...
            gamma, sub, num_subjects);
        
        % Initialize matrices to store DTW distances for the current subject.
        % These matrices are of size [num_components x num_components] where
        % only the lower triangular part (unique component pairs) is filled.
        sub_tr_dtw     = zeros(num_timepoints, num_components, num_components);
        sub_dtw_fii     = zeros(num_timepoints-win_size, num_components, num_components);
        
        % Loop over all unique component pairs (comp1 > comp2)
        for comp1 = 2:num_components
            for comp2 = 1:(comp1 - 1)
                % Extract the time courses for the two components and z-score them.
                x = squeeze(subTcs(sub, :, comp1));
                y = squeeze(subTcs(sub, :, comp2));
                x = zscore(x);
                y = zscore(y);
                
                % Compute time-resolved DTW using the custom DTW framework.
                [~, ~, tr_dtw, ~] = DTWFramework(x, y, gamma, win_size);
                
                % Store the computed DTW distances in the lower-triangular matrix.
                sub_tr_dtw(:, comp1, comp2) = tr_dtw;

                sub_dtw_cumm = zeros(num_cumm, 1);
                
                for t = t_a:num_timepoints
                    sub_dtw_cumm(t - t_a + 1) = mean(tr_dtw(1:t));
                end
                
                smoothed_cumm = gauss_smooth(sub_dtw_cumm, FWHM);
                sub_dtw_fii(:, comp1, comp2) = diff(smoothed_cumm);
            end
        end
        
        % Vectorize the lower-triangular matrices and store them in the output arrays.
        full_tr_dtw(sub, :, :) = icatb_mat2vec(sub_tr_dtw);
        full_dtw_fii(sub, :, :) = icatb_mat2vec(sub_dtw_fii);

    end
end
