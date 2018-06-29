% Comparison of G8/G28 Hall MHD and MHD-EPIC reconnection efficiencies.
%
% Require pre-processed CPCP data.
%
% Hongyang Zhou, hyzhou@umich.edu 02/26/2018

clear; clc
%% Read data and preprocess
%G8 = load('CPCP_G8.mat');
G8 = load('~/Ganymede/newPIC/CPCPdata_theta51/CPCP_G8_51.mat');
G28 = load('CPCP_G28.mat');
G8Hall = load('CPCP_G8Hall.mat');
G28Hall = load('CPCP_G28Hall.mat');
%G2 = load('CPCP_G2.mat');

% Remove the duplicate time(s)
G28.time(601) = [];
G28.CPCPt(601) = [];
G28.Potential_bk(601) = [];

rateG8 = G8.CPCPt ./ G8.Potential_bk;
rateG28 = G28.CPCPt ./ G28.Potential_bk;
rateG8Hall = G8Hall.CPCPt ./ G8Hall.Potential_bk;
rateG28Hall = G28Hall.CPCPt ./ G28Hall.Potential_bk;

% Ignore the transition time (100s) for FFT calculation
TimeTransition = 100; %[s]

rateMeanG8  = mean(rateG8(1:TimeTransition:end));
rateMeanG28 = mean(rateG28(1:TimeTransition:end));
rateMeanG8Hall  = mean(rateG8Hall(1:TimeTransition:end));
rateMeanG28Hall = mean(rateG28Hall(1:TimeTransition:end));

% FFT
Fs = 1; % Sample rate is 1 per second
NFFT = length(rateG8) - TimeTransition; % Number of FFT points
freq = (0 : 1/NFFT : 1/2-1/NFFT)*Fs; % Frequency vector 
period = 1./freq;

fftG8 = fft(rateG8(1+TimeTransition:end),NFFT);
fftG8(1) = 0;   % Remove the DC component for better visualization
powerG8 = abs(fftG8(1:floor(NFFT/2))).^2;

fftG28 = fft(rateG28(1+TimeTransition:end),NFFT);
fftG28(1) = 0;  % Remove the DC component for better visualization
powerG28 = abs(fftG28(1:floor(NFFT/2))).^2;

fftG8Hall = fft(rateG8Hall(1+TimeTransition:end),NFFT);
fftG8Hall(1) = 0; % Remove the DC component for better visualization
powermix1 = abs(fftG8Hall(1:floor(NFFT/2))).^2;

fftG28Hall = fft(rateG28Hall(1+TimeTransition:end),NFFT);
fftG28Hall(1) = 0; % Remove the DC component for better visualization
powermix2 = abs(fftG28Hall(1:floor(NFFT/2))).^2;

%% Reconnection rate
addpath('~/Ganymede/scripts/PlotPub/lib');

% G8
figure;
plot(G8.time,rateG8,G8Hall.time,rateG8Hall)
opt = [];
opt.XLabel = 'Simulations time [s]';
opt.YLabel = 'Reconnection Efficiency';
opt.Legend = {'G8','G8Hall'};
opt.BoxDim = [10, 5]; %[width, height]

setPlotProp(opt)

% G28
figure;
plot(G28.time,rateG28,G28Hall.time,rateG28Hall)
opt = [];
opt.XLabel = 'Simulations time [s]';
opt.YLabel = 'Reconnection Efficiency';
opt.Legend = {'G28','G28Hall'};
opt.BoxDim = [10, 5]; %[width, height]

setPlotProp(opt)
