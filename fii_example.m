%% Example Script: Computing Functional Inertia Index (FII) for a Single Pair
clc; clear;
addpath(genpath('./'));

%% Parameters
Ts       = 2;                 % Sampling time (same as TR)
Fs       = 1/Ts; 
band     = [0.01 0.15];       % BOLD band (Hz)
cutoff   = [0.6 1.5];         % Filter stability parameters
num_timepoints = 1000;

gamma   = 1.5;                % DTW exponent (from test–retest optimization)
win_size = 45;                % Choose an anchor window (e.g. 45) in samples
FWHM    = 5;                  % Gaussian smoothing before derivative

%% Generate two example band-limited random signals
x = bandpass_filtering(randn(1, num_timepoints), Ts, band, cutoff);
y = bandpass_filtering(randn(1, num_timepoints), Ts, band, cutoff);

% Normalize as done in the main pipeline
x = zscore(x);
y = zscore(y);

%% Compute DTW-FII
[dtw_fii, tr_dtw] = compute_dtw_fii(x, y, gamma, FWHM, win_size);

% Optional: flip sign so positive = stabilizing (as used in the manuscript)
dtw_fii = -dtw_fii;

%% Plot results
figure;

subplot(2,1,1)
t_fii = win_size+1:num_timepoints;   % FII starts after anchor window
plot(t_fii, dtw_fii, 'LineWidth', 1.5)
title('Functional Inertia Index (FII)')
xlabel('Time (samples)')
ylabel('FII')
grid on

subplot(2,1,2)
plot(tr_dtw, 'LineWidth', 1)
title('Time-Resolved DTW (Amplitude Discrepancy)')
xlabel('Time (samples)')
ylabel('tr-DTW')
grid on