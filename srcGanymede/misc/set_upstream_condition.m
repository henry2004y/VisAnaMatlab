% Ganymede time-accurate simulation upstream boundary conditions
%
% Hongyang Zhou, 05/24/2019

clear; clc; close all
%%
nPoints = 100;
t = linspace(0,10,nPoints);
RotationPeriod = 10*3600; % [s]
omega = 2*pi/RotationPeriod; %
phase = linspace(0,360,nPoints); % West Longitude
offset = [-115,0,5,-115+90+180,0];

% Baseline magnetic field
B0 = [0, 70, 0];
% Perturbation
Mag = [83, 0, 10];

%% Compare with KK97 analytical model
% Figure 2, Jia et.al [2008]
Br = B0(1) + Mag(1).*sind(phase+offset(1));
Bt = B0(2) + Mag(2).*sind(phase+offset(2));
Bp = B0(3) + Mag(3).*sind(phase+offset(3));

figure
subplot(311)
plot(phase,Br,'LineWidth',1.2)
yline(0)
xlim([0 360])
ylabel('Br [nT]')
set(gca,'XMinorTick','on','LineWidth',1.2,'FontSize',14)
subplot(312)
plot(phase,Bt,'LineWidth',1.2)
xlim([0 360])
ylabel('B\theta [nT]')
set(gca,'XMinorTick','on','LineWidth',1.2,'FontSize',14)
subplot(313)
plot(phase,Bp,'LineWidth',1.2)
xlim([0 360])
ylabel('B\phi [nT]')
xlabel('West Longitude')
set(gca,'XMinorTick','on','LineWidth',1.2,'FontSize',14)


%% Transform into GPhiO coordinate system

Bx = Bp;
By = -Br;
Bz = -Bt;

% Density
n = 3 + 1.0.*sind(phase+offset(4)); %[amu/cc]

% Pressure, assuming constant temperature
kb = 1.38e-23;
T = 5.7367e7;    %[K]
P = n*kb*T*1e15; %[nPa]


%%
figure('Position',[89 29 749 771])
subplot(511)
plot(phase,Bx,'LineWidth',1.2)
xlim([0 360])
ylabel('Bx [nT]')
set(gca,'XMinorTick','on','LineWidth',1.2,'FontSize',14)
subplot(512)
plot(phase,By,'LineWidth',1.2)
xlim([0 360])
ylabel('By [nT]')
set(gca,'XMinorTick','on','LineWidth',1.2,'FontSize',14)
subplot(513)
plot(phase,Bz,'LineWidth',1.2)
xlim([0 360])
ylabel('Bz [nT]')
set(gca,'XMinorTick','on','LineWidth',1.2,'FontSize',14)
subplot(514)
plot(phase,n,'LineWidth',1.2)
xlim([0 360])
ylabel('n [amu/cc]')
set(gca,'XMinorTick','on','LineWidth',1.2,'FontSize',14)
subplot(515)
plot(phase,P,'LineWidth',1.2)
xlim([0 360])
ylabel('P [nPa]')
set(gca,'XMinorTick','on','LineWidth',1.2,'FontSize',14)
xlabel('West Longitude')

%% Generate input file

% File format:
% Year Month Day Hour Min Sec Msec Bx[nT] By[nT] Bz[nT] 
% Vx[km/s] Vy[km/s] Vz[km/s] N[cm^(-3)] T[Kelvin]

fileID = fopen('UpstreamInput.txt','w');

fmt = ['%4d %4d %4d %4d %4d %4d %4d %10.5f %10.5f %10.5f ' ...
   '%10.3f %10.3f %10.3f %10.3f %14.2f\n'];

yr = Inf(1,nPoints);
month = Inf(1,nPoints);
day = Inf(1,nPoints);
hr = Inf(1,nPoints);
mn = Inf(1,nPoints);
sec = Inf(1,nPoints);
msec = zeros(1,nPoints); % As an approximation, don't count msec
Vx = ones(1,nPoints); % GSM
Vy = zeros(1,nPoints);
Vz = zeros(1,nPoints);
T  = 5.7367e7.*ones(1,nPoints);

Vx(:) = 140; %[km/s]

tStart = datetime([2000  1  1  0  0  0]);
tEnd   = tStart + seconds(RotationPeriod);
time   = linspace(tStart,tEnd,nPoints);

Input = [time.Year; time.Month; time.Day; time.Hour; ...
   time.Minute; round(time.Second); msec; ...
   Bx; By; Bz; Vx; Vy; Vz; n; T];


fprintf(fileID,fmt,Input);
     
fclose(fileID);
     