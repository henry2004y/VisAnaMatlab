% Script for SWMF Galileo Flyby Simulation Comparison 
% Purpose: Plot time-accurate magnetic field comparisons.
%
% Prerequisite: Folder FileIO_SWMF must be included into Matlab path
%
% This script is able to handle any npict number of snapshots. If the 
% simulation time interval is smaller than observation (which is usually
% the case), than the magnetic field along the trajectoryy before 
% time_start is obtained from first single snapshot ipict=1, and B after
% time_end is obtained from the last single snapshot ipict=npict. In this 
% way, I avoid the jump that Gabor had due to periodic interpolation.
%
% flyby and time_start need to be modified when switching between flybys.
%
% Hongyang Zhou, hyzhou@umich.edu 01/07/2020

%% Parameters
flyby =  8;   % [1,2,7,8,28,29]
DoPlot = 1;  % Plot output
DoSave = 0;  % Save norm2 number
firstpict = 1; % first snapshot to pick
lastpict  = 1; % last snapshot to pick
fileGathered = false; % Input data in 1 file or multiple files

%% Read observation data
flybyfile = strcat('Galileo_G',int2str(flyby),'_flyby_MAG.dat');
f = fullfile('~/Ganymede/GalileoData/galileomagdata',flybyfile);
[~,data] = read_log_data(f);

time = datetime(data(:,1:6));
xyz  = data(:,7:9);
Bobs = data(:,10:12);
BobsStrength = sqrt(Bobs(:,1).^2+Bobs(:,2).^2+Bobs(:,3).^2);

%% Read/Plot simulation data 
switch flyby
   case 8
      %filename='~/Ganymede/MOP2018/runG8_PIC_1200s/box_B_Galileo.outs';
      %filename='~/Ganymede/Hall_AMR3/GM/box_Hall_B_1200.outs';
      %filename='~/Ganymede/Hall_AMR2/box_900s.outs';
      filename='~/Ganymede/PIC_frontera/GM/box_PIC_B_1200.outs';
      % Select the starting time for synthetic satellite
      time_start = datetime(1997,5,7,15,47,30);
      %time_start = datetime(1997,5,7,15,45,40); % G8 PIC
      %time_start = datetime(1997,5,7,15,48,7); % G8
   case 28
      filename='~/Ganymede/MOP2018/runG28_PIC_1200s/GM/box_B_Galileo.outs';
      % Select the starting time for synthetic satellite
      %time_start = datetime(2000,5,20,9,48,0); % G28
      %time_start = datetime(2000,5,20,9,52,35); % G28
      %time_start = datetime(2000,5,20,9,51,40); % G28
      time_start = datetime(2000,5,20,10,0,12);
end
[filehead,~,filelist] = read_data(filename,'verbose',false);
nx = filehead.nx; nw = filehead.nw;
npict = filelist.npictinfiles;

if firstpict>npict
   error('firstpict out of range!')
end

npict = npict + 1 - firstpict;

%Tobs = datenum(time); dt = 1/(24*3600);
%Timerange_obs = Tobs(end)-Tobs(1);
% num_sim = floor(Timerange_obs/dt);

kStart = find( time_start-time<0,1 );
%kStart = find( time_start-time>0,1,'last' );

time_end   = time_start + (npict-1)/(24*3600);
kEnd = find( time_end-time<0,1 );

num_sim = kStart + (numel(time) - kEnd + 1) + npict - 2;
Bsim = Inf(num_sim,3);
timesim = NaT(num_sim,1);
timesim(1:kStart) = time(1:kStart);
% This is not that accurate, but just an approximation
timesim(kStart:kStart+npict-1) = linspace(time_start,time_end,npict);
timesim(kStart+npict-1:end) = time(kEnd:end);

% Extract the magnetic field along the trajectory before the starting time
% in one snapshot
[~,data] = read_data(filename,'npict',firstpict,'verbose',false);

% Interpolate simulation data to observation data
data = data.file1;

x = data.x(:,:,:,1);
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

w = data.w;
bx = data.w(:,:,:,1);
by = data.w(:,:,:,2);
bz = data.w(:,:,:,3);

% From ndgrid to meshgrid format
% x  = permute(x,[2 1 3]);
% y  = permute(y,[2 1 3]);
% z  = permute(z,[2 1 3]);
% bx = permute(bx,[2 1 3]);
% by = permute(by,[2 1 3]);
% bz = permute(bz,[2 1 3]);

% Bsim(1:kStart,1) = interp3(x,y,z,bx,xyz(1:kStart,1),xyz(1:kStart,2),xyz(1:kStart,3));
% Bsim(1:kStart,2) = interp3(x,y,z,by,xyz(1:kStart,1),xyz(1:kStart,2),xyz(1:kStart,3));
% Bsim(1:kStart,3) = interp3(x,y,z,bz,xyz(1:kStart,1),xyz(1:kStart,2),xyz(1:kStart,3));

F1 = griddedInterpolant(x,y,z,bx);
F1.Method = 'linear';
F2 = griddedInterpolant(x,y,z,by);
F2.Method = 'linear';
F3 = griddedInterpolant(x,y,z,bz);
F3.Method = 'linear';

Bsim(1:kStart,1) = F1(xyz(1:kStart,1),xyz(1:kStart,2),xyz(1:kStart,3));
Bsim(1:kStart,2) = F2(xyz(1:kStart,1),xyz(1:kStart,2),xyz(1:kStart,3));
Bsim(1:kStart,3) = F3(xyz(1:kStart,1),xyz(1:kStart,2),xyz(1:kStart,3));


% From ndgrid to meshgrid format
x  = permute(x,[2 1 3]);
y  = permute(y,[2 1 3]);
z  = permute(z,[2 1 3]);
bx = permute(bx,[2 1 3]);
by = permute(by,[2 1 3]);
bz = permute(bz,[2 1 3]);
   
% Extract the magnetic field along the trajectory from multiple snapshot
for ipict=1:npict-1
   fprintf('ipict=%d\n',ipict);
   [~,data] = read_data(filename,'npict',firstpict+ipict,'verbose',false);
   % Interpolate simulation data to observation data
   data = data.file1;
   
   x = data.x(:,:,:,1);
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);
   
   w = data.w;
   bx = data.w(:,:,:,1);
   by = data.w(:,:,:,2);
   bz = data.w(:,:,:,3);
   
   % From ndgrid to meshgrid format
   x  = permute(x,[2 1 3]);
   y  = permute(y,[2 1 3]);
   z  = permute(z,[2 1 3]);
   bx = permute(bx,[2 1 3]);
   by = permute(by,[2 1 3]);
   bz = permute(bz,[2 1 3]);
   
   [~,k] = min( abs(timesim(kStart+ipict) - time) );
   
   Bsim(kStart+ipict,1) = interp3(x,y,z,bx,xyz(k,1),xyz(k,2),xyz(k,3));
   Bsim(kStart+ipict,2) = interp3(x,y,z,by,xyz(k,1),xyz(k,2),xyz(k,3));
   Bsim(kStart+ipict,3) = interp3(x,y,z,bz,xyz(k,1),xyz(k,2),xyz(k,3));
   
end

% Extract the magnetic field along the trajectory after the ending time
% in one snapshot 
% Bsim(kStart+npict:end,1) = interp3(x,y,z,bx,xyz(kEnd+1:end,1),xyz(kEnd+1:end,2),xyz(kEnd+1:end,3));
% Bsim(kStart+npict:end,2) = interp3(x,y,z,by,xyz(kEnd+1:end,1),xyz(kEnd+1:end,2),xyz(kEnd+1:end,3));
% Bsim(kStart+npict:end,3) = interp3(x,y,z,bz,xyz(kEnd+1:end,1),xyz(kEnd+1:end,2),xyz(kEnd+1:end,3));

% From ndgrid to meshgrid format
x  = permute(x,[2 1 3]);
y  = permute(y,[2 1 3]);
z  = permute(z,[2 1 3]);
bx = permute(bx,[2 1 3]);
by = permute(by,[2 1 3]);
bz = permute(bz,[2 1 3]);

   F1 = griddedInterpolant(x,y,z,bx);
   F1.Method = 'linear';
   F2 = griddedInterpolant(x,y,z,by);
   F2.Method = 'linear';
   F3 = griddedInterpolant(x,y,z,bz);
   F3.Method = 'linear';
   
   Bsim(kStart+npict:end,1) = F1(xyz(kEnd+1:end,1),xyz(kEnd+1:end,2),xyz(kEnd+1:end,3));
   Bsim(kStart+npict:end,2) = F2(xyz(kEnd+1:end,1),xyz(kEnd+1:end,2),xyz(kEnd+1:end,3));
   Bsim(kStart+npict:end,3) = F3(xyz(kEnd+1:end,1),xyz(kEnd+1:end,2),xyz(kEnd+1:end,3)); 

BsimStrength = sqrt(Bsim(:,1).^2+Bsim(:,2).^2+Bsim(:,3).^2);   
   
%% Visualization

% I can try plot(Bsim) to put all three components into one plot

if DoPlot
   figure('Position', [100, 100, 1000, 700]);
   tStart_ = 2000;
   tEnd_   = 5000;
   h = subplot(411); LW1 = 2.5; LW2 = 1.5; FS=16; Height = 0.2;
   plot(time, Bobs(:,1),'k',timesim,Bsim(:,1),'r','LineWidth',LW1); %Bx
   h.XTickLabel = []; 
   h.Position(4) = Height;
   ylabel('Bx [nT]');
   legend(h,{'Obs','Sim'})
   xlim([time(tStart_) time(tEnd_)]);
   ylim([min(Bobs(:,1))-50 max(Bobs(:,1))+50 ]);
   %title({'(a)',strcat('Galileo G',int2str(flyby),' Flyby Magnetic field')})
   title(strcat('Galileo G',int2str(flyby),' Flyby Magnetic field'))
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   
   %hold on
   %xline(time(4000),'--r');
   
   h = subplot(412);
   plot(time, Bobs(:,2),'k',timesim,Bsim(:,2),'r','LineWidth',LW1); %By
   h.XTickLabel = [];
   h.Position(4) = Height;   
   ylabel('By [nT]');
   xlim([time(tStart_) time(tEnd_)]);
   ylim([min(Bobs(:,2))-50 max(Bobs(:,2))+50 ]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(413);
   plot(time, Bobs(:,3),'k',timesim, Bsim(:,3),'r','LineWidth',LW1); %Bz
   h.XTickLabel = [];
   h.Position(4) = Height;   
   ylabel('Bz [nT]');
   xlim([time(tStart_) time(tEnd_)]);
   ylim([min(Bobs(:,3))-50 max(Bobs(:,3))+50 ]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(414);   
   plot(time,BobsStrength,'k',timesim,BsimStrength,'r','LineWidth',LW1); %|B|
   h.Position(4) = Height;   
   ylabel('B [nT]');
   xlim([time(tStart_) time(tEnd_)]);
   ylim([min(BobsStrength)-50 max(BobsStrength)+50 ]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
end

%% Save simulation line data
if DoSave
   save('B_G8_Hall_AMR2.mat','Bobs','Bsim','time','timesim')
end

%% For the App

%tmt = timetable(time, Bobs, BobsStrength);
%tmt.Properties.VariableNames = ["Bobs" "B"];

%clear data

%%
fs = 3; % [Hz]
window = 128;
nOverlap = window / 2;

figure
subplot(1,2,1)
spectrogram(BobsStrength(1750:6295), window, nOverlap, [], fs, 'yaxis')
%ax = gca;
%ax.YScale = 'log';
ylim([0,0.5])

fs = 1; % [Hz]
window = 40;
nOverlap = 0;

subplot(1,2,2)
spectrogram(BsimStrength(1894:3408), window, nOverlap, [], fs, 'yaxis')
%ylim([0,0.5])
yt = get(gca, 'YTick');
set(gca, 'YTick', 0:0.05:0.5, 'YTickLabel', (0:0.05:0.5)*1e3);

%%

% Check https://www.mathworks.com/help/matlab/ref/tiledlayout.html
% for arranging the subplots!
% tile only works after version 2019b!

t = tiledlayout(2,2,'TileSpacing','None');

%figure('Position', [100, 100, 1000, 700]);
%ax(1) = subplot(221);
nexttile(1)
range_obs_ = 1600:6200;
plot(time(range_obs_),BobsStrength(range_obs_),'k','LineWidth',LW1); %|B|
ylim([0,200])
title('Observation B')   

%ax(2) = subplot(222);
nexttile(2)
range_sim_ = 1600:3400;
plot(timesim(range_sim_),BsimStrength(range_sim_),'r','LineWidth',LW1); %|B|
ylim([0,200])
title('Simulation B')

%ax(3) = subplot(223);
nexttile(3)
fs = 3; % [Hz]
window = 128;
nOverlap = window / 2;
[S,F,T] = spectrogram(BobsStrength(1750:6295), window, nOverlap, [], fs, 'yaxis');
surf(T/60,F,abs(S), 'edgecolor', 'none')
colorbar
ylim([0,0.5])
yticks(0:0.05:0.5)
%shading interp
caxis([1, 11621])
set(gca,'ColorScale','log')
view([0 90])
axis tight
xlabel('Time [s]')
ylabel('Frequency [Hz]')
%set(gca,'YScale','log')

%ax(4) = subplot(224);
nexttile(4)
fs = 1; % [Hz]
window = 40;
nOverlap = window / 2;
[S,F,T] = spectrogram(BsimStrength(1894:3408), window, nOverlap, [], fs, 'yaxis');
surf(T/60,F,abs(S), 'edgecolor', 'none')
colorbar
caxis([1, 11621])
colormap bone
%shading interp
set(gca,'ColorScale','log')
view([0 90])
axis tight
xlabel('Time')
ylabel('Frequency')
%set(gca,'YScale','log')
