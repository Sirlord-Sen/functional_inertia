function [idx, c, sumd, D] = calculate_kmeans(estimation_2D, cluster_num, distance, rep, maxIter, useParallel, doDisplay)
%CALCULATE_KMEANS  K-means clustering with optional parameters
%
% Usage:
%   [idx,c,sumd,D] = calculate_kmeans(X, k)
%   [idx,c,sumd,D] = calculate_kmeans(X, k, distance)
%   [idx,c,sumd,D] = calculate_kmeans(X, k, distance, rep)
%   [idx,c,sumd,D] = calculate_kmeans(X, k, distance, rep, maxIter)
%   [idx,c,sumd,D] = calculate_kmeans(X, k, distance, rep, maxIter, useParallel)
%   [idx,c,sumd,D] = calculate_kmeans(X, k, distance, rep, maxIter, useParallel, doDisplay)
%
% Inputs:
%   estimation_2D : (n×m) data matrix to cluster
%   cluster_num   : number of clusters
%   distance      : distance metric (string). Default: 'sqeuclidean'
%   rep           : number of replicates. Default: 20
%   maxIter       : maximum iterations. Default: 10000
%   useParallel   : flag (0 or 1) to enable parallel and substreams. Default: 1
%   doDisplay     : flag (0 or 1) to show kmeans output. Default: 1 (display)
%
% Outputs:
%   idx   : cluster indices for each row of X
%   c     : cluster centroids
%   sumd  : within-cluster sums of point-to-centroid distances
%   D     : distances from each point to every centroid

    % Set defaults for optional inputs
    if nargin < 3 || isempty(distance)
        distance = 'sqeuclidean';
    end
    if nargin < 4 || isempty(rep)
        rep = 20;
    end
    if nargin < 5 || isempty(maxIter)
        maxIter = 10000;
    end
    if nargin < 6 || isempty(useParallel)
        useParallel = 0;
    end
    if nargin < 7 || isempty(doDisplay)
        doDisplay = 0;
    end

    % Configure parallel options
    if useParallel
        stream = RandStream('mlfg6331_64');
        options = statset('UseParallel',1, 'UseSubstreams',1, 'Streams',stream);
    else
        options = statset('UseParallel',0, 'UseSubstreams',0);
    end

    % Set display option based on doDisplay flag
    if doDisplay
        displayMode = 'final';  % or 'iter' if you want full progress
    else
        displayMode = 'off';
    end

    % Run kmeans
    [idx, c, sumd, D] = kmeans(estimation_2D, cluster_num, ...
                               'Distance', distance, ...
                               'Replicates', rep, ...
                               'MaxIter', maxIter, ...
                               'Options', options, ...
                               'Display', displayMode);
end
