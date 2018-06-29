% Plotting the CPCP-time relation and frequency analysis.
%
% Could be used after the calculation of CPCP for all four runs.
%
% Requirement:
% PlotPub must be added to path.
%
% Hongyang Zhou, hyzhou@umich.edu 01/11/2018

clear;clc; %close all; 
%% Read data and preprocess
% G8 = load('~/Ganymede/newPIC/CPCPdata_theta50/CPCP_G8.mat');
% G28 = load('~/Ganymede/newPIC/CPCPdata_theta50/CPCP_G28.mat');
mix1 = load('~/Ganymede/newPIC/CPCPdata_theta50/CPCP_mix1.mat');
mix2 = load('~/Ganymede/newPIC/CPCPdata_theta50/CPCP_mix2.mat');

% Remove the duplicate time(s)
% G28.time(601) = [];
% G28.CPCPt(601) = [];
% G28.Potential_bk(601) = [];

mix1.time(601) = [];
mix1.CPCPt(601) = [];
mix1.Potential_bk(601) = [];

mix2.time(601) = [];
mix2.CPCPt(601) = [];
mix2.Potential_bk(601) = [];


G8 = load('~/Ganymede/newPIC/CPCPdata_theta51/CPCP_G8_51.mat');
G8.time(301) = [];
G8.CPCPt(301) = [];
G8.Potential_bk(301) = [];
G8.time(601) = [];
G8.CPCPt(601) = [];
G8.Potential_bk(601) = [];
G8.time(901) = [];
G8.CPCPt(901) = [];
G8.Potential_bk(901) = [];
G28 = load('~/Ganymede/newPIC/CPCPdata_theta51/CPCP_G28_51.mat');

% G8 = load('CPCP_G8.mat');
% G28 = load('CPCP_G28.mat');
% mix1 = load('CPCP_mix1.mat');
% mix2 = load('CPCP_mix2.mat');



% Ignore the transition time (100s) for FFT calculation
TimeTransition = 100; %[s]

%% Test
% use reconnection to do fft instead of CPCP itself

rateG8 = G8.CPCPt ./ G8.Potential_bk;
rateG28 = G28.CPCPt ./ G28.Potential_bk;
ratemix1 = mix1.CPCPt ./ mix1.Potential_bk;
ratemix2 = mix2.CPCPt ./ mix2.Potential_bk;

rateMeanG8  = mean(rateG8(1:TimeTransition:end));
rateMeanG28 = mean(rateG28(1:TimeTransition:end));
rateMeanMix1= mean(ratemix1(1:TimeTransition:end));
rateMeanMix2= mean(ratemix2(1:TimeTransition:end));

% Smooth the original data
% rateG8Smooth   = smoothdata(rateG8);
% rateG28Smooth  = smoothdata(rateG28);
% ratemix1Smooth = smoothdata(ratemix1);
% ratemix2Smooth = smoothdata(ratemix2);

%% FFT
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

fftmix1 = fft(ratemix1(1+TimeTransition:end),NFFT);
fftmix1(1) = 0; % Remove the DC component for better visualization
powermix1 = abs(fftmix1(1:floor(NFFT/2))).^2;

fftmix2 = fft(ratemix2(1+TimeTransition:end),NFFT);
fftmix2(1) = 0; % Remove the DC component for better visualization
powermix2 = abs(fftmix2(1:floor(NFFT/2))).^2;

% FFT of smoothed data (This is actually wrong!!!)
% fftG8Smooth   = fft(rateG8Smooth(1+TimeTransition:end),NFFT);
% fftG28Smooth  = fft(rateG28Smooth(1+TimeTransition:end),NFFT);
% fftmix1Smooth = fft(ratemix1Smooth(1+TimeTransition:end),NFFT);
% fftmix2Smooth = fft(ratemix2Smooth(1+TimeTransition:end),NFFT);
% 
% powerG8Smooth   = abs(fftG8Smooth(1:floor(NFFT/2))).^2;
% powerG28Smooth  = abs(fftG28Smooth(1:floor(NFFT/2))).^2;
% powermix1Smooth = abs(fftmix1Smooth(1:floor(NFFT/2))).^2;
% powermix2Smooth = abs(fftmix2Smooth(1:floor(NFFT/2))).^2;

%% Visualization

%% Global reconnection rate
figure;
plot(G8.time,G8.CPCPt,G28.time,G28.CPCPt,mix1.time,mix1.CPCPt,...
   mix2.time,mix2.CPCPt)
legend({'G8','G28','mix1','mix2'})
xlabel('Simulations time [s]');
ylabel('Total reconnection rate');
set(gca,'FontSize',16,'LineWidth',1.2);

%% Global reconnection efficiency
figure;
plot(G8.time,rateG8,G28.time,rateG28,mix1.time,ratemix1,...
   mix2.time,ratemix2)
legend({'G8','G28','mix1','mix2'})
xlabel('Simulations time [s]');
ylabel('Reconnection Efficiency');
set(gca,'FontSize',16,'LineWidth',1.2);

%%
% better quality for the paper
figure;
addpath('~/Ganymede/scripts/PlotPub/lib');
plot(G8.time,rateG8,G28.time,rateG28,mix1.time,ratemix1,...
   mix2.time,ratemix2)
opt = [];
opt.XLabel = 'Simulations time [s]';
opt.YLabel = 'Reconnection Efficiency';
opt.Legend = {'G8','G28','mix1','mix2'};
opt.BoxDim = [10, 5]; %[width, height]

setPlotProp(opt)

%%
% Frequency domain amplitude
% figure('Position', [100, 100, 1049, 895]);
% LW1 = 3; LW2 = 1.5; FS=16;
% subplot(411)
% plot(freq,abs(fftG8(1:floor(NFFT/2))))
% set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
% subplot(412)
% plot(freq,abs(fftG28(1:floor(NFFT/2))))
% set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
% subplot(413)
% plot(freq,abs(fftmix1(1:floor(NFFT/2))))
% set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
% subplot(414)
% plot(freq,abs(fftmix2(1:floor(NFFT/2))))
% set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
% xlabel('Frequency (cycles/second)');
% ylabel('Magnitude')

%% Frequency domain log-log plot
figure('Position', [100, 100, 1049, 895]);
LW1 = 3; LW2 = 1.5; FS=16;
subplot(411)
loglog(freq,powerG8)
xlim([1e-3 1]); ylim([1e-5 1e5])
set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
subplot(412)
loglog(freq,powerG28)
xlim([1e-3 1]); ylim([1e-5 1e5])
set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
subplot(413)
loglog(freq,powermix1)
xlim([1e-3 1]); ylim([1e-5 1e5])
set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
subplot(414)
loglog(freq,powermix2)
xlim([1e-3 1]); ylim([1e-5 1e5])
set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
xlabel('Frequency (cycles/second)');
ylabel('Magnitude')


%% Frequency domain log-log plot for smoothed data in one plot
% figure('Position', [100, 100, 1000, 850]);
% LW1 = 3; LW2 = 1.5; FS=16;
% loglog(freq,powerG8Smooth,freq,powerG28Smooth,freq,powermix1Smooth,...
%    freq,powermix2Smooth)
% 
% legend({'G8','G28','mix1','mix2'})
% xlabel('Frequency (cycles/second)');
% ylabel('Power')
% set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);

% Better quality for the paper
figure;
loglog(freq,powerG8,freq,powerG28,freq,powermix1,freq,powermix2)
opt = [];
opt.XLabel = 'Frequency (cycles/second)';
opt.YLabel = 'Power';
opt.Legend = {'G8','G28','mix1','mix2'};
opt.BoxDim = [10, 4]; %[width, height]

setPlotProp(opt)

% figure;
% loglog(freq,powerG8Smooth,freq,powerG28Smooth,freq,powermix1Smooth,...
%    freq,powermix2Smooth)
% opt = [];
% opt.XLabel = 'Frequency (cycles/second)';
% opt.YLabel = 'Power';
% opt.Legend = {'G8','G28','mix1','mix2'};
% opt.BoxDim = [10, 4]; %[width, height]
% 
% setPlotProp(opt)


%% Switch to period power plot
% Usually people use [dbw], which a log scale power  
figure('Position', [100, 100, 1049, 895]);
LW1 = 3; LW2 = 1.5; FS=16;
subplot(411)
%plot(period,powerG8)
semilogy(period,powerG8)
ylabel('G8')
set(gca,'FontSize',FS,'xMinorTick','on','yMinorTick','on',...
   'LineWidth',LW2);
xlim([0 400])
subplot(412)
%plot(period,powerG28)
semilogy(period,powerG28)
ylabel('G28')
set(gca,'FontSize',FS,'xMinorTick','on','yMinorTick','on',...
   'LineWidth',LW2);
xlim([0 400])
subplot(413)
%plot(period,powermix1)
semilogy(period,powermix1)
ylabel('mix1')
set(gca,'FontSize',FS,'xMinorTick','on','yMinorTick','on',...
   'LineWidth',LW2);
xlim([0 400])
subplot(414)
%plot(period,powermix2)
semilogy(period,powermix2)
ylabel('mix2')
set(gca,'FontSize',FS,'xMinorTick','on','yMinorTick','on',...
   'LineWidth',LW2);
xlim([0 400])
xlabel('Period (second/cycle)');


%% FIR Filter
% Create ?sgolayfilt? Filtered FFT
sgfG8 = sgolayfilt(powerG8, 1, 61);
sgfG28 = sgolayfilt(powerG28, 1, 41);
sgfmix1 = sgolayfilt(powermix1, 1, 61);
sgfmix2 = sgolayfilt(powermix2, 1, 41);

% figure
% loglog(freq,powermix2); hold on
% loglog(freq,sgfmix2,'-r','LineWidth',1.2); hold off

figure;
loglog(freq,sgfG8,freq,sgfG28,freq,sgfmix1,freq,sgfmix2)
opt = [];
opt.XLabel = 'Frequency (cycles/second)';
opt.YLabel = 'Power';
opt.Legend = {'G8','G28','mix1','mix2'};
opt.BoxDim = [10, 4]; %[width, height]

setPlotProp(opt)


% Power series fit
k = find(sgfG8<4.7);
fG8 = fit(freq(k)',sgfG8(k),'power1');
k = find(sgfG28<3.9);
fG28 = fit(freq(k)',sgfG28(k),'power1');
k = find(sgfmix1<2.4);
fmix1 = fit(freq(k)',sgfmix1(k),'power1');
k = find(sgfmix2<1.5);
fmix2 = fit(freq(k)',sgfmix2(k),'power1');
figure;
loglog(freq,sgfG8,freq,sgfG28,freq,sgfmix1,freq,sgfmix2); hold on
plot(fG8)
plot(fG28)
plot(fmix1)
plot(fmix2)
