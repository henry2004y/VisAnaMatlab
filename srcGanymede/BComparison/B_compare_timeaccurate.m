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
% Hongyang Zhou, hyzhou@umich.edu 01/03/2018

%% Parameters
flyby = 28;   % [1,2,7,8,28,29]
DoPlot = 1;  % Plot output
DoSave = 0;  % Save norm2 number
firstpict = 61; % first snapshot to pick
lastpict  = 1; % last snapshot to pick

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
      filename='~/Ganymede/newPIC/run_G8_newPIC/box_B_G8_1200s.outs';
      % Select the starting time for synthetic satellite
      time_start = datetime(1997,5,7,15,45,2); % G8
   case 28
      filename='~/Ganymede/newPIC/run_G28_newPIC/box_B_G28_1200s.outs';
      % Select the starting time for synthetic satellite
      %time_start = datetime(2000,5,20,9,48,0); % G28
      %time_start = datetime(2000,5,20,9,52,35); % G28
      time_start = datetime(2000,5,20,9,51,40); % G28
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
   h = subplot(411); LW1 = 2.5; LW2 = 1.5; FS=16; Height = 0.2;
   plot(time, Bobs(:,1),'k',timesim,Bsim(:,1),'r','LineWidth',LW1); %Bx
   h.XTickLabel = []; 
   h.Position(4) = Height;
   ylabel('Bx [nT]');
   legend(h,{'Obs','Sim'})
   xlim([min(time) max(time)]);
   ylim([min(Bobs(:,1))-50 max(Bobs(:,1))+50 ]);
   %title({'(a)',strcat('Galileo G',int2str(flyby),' Flyby Magnetic field')})
   title(strcat('Galileo G',int2str(flyby),' Flyby Magnetic field'))
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(412);
   plot(time, Bobs(:,2),'k',timesim,Bsim(:,2),'r','LineWidth',LW1); %By
   h.XTickLabel = [];
   h.Position(4) = Height;   
   ylabel('By [nT]');
   xlim([min(time) max(time)]);
   ylim([min(Bobs(:,2))-50 max(Bobs(:,2))+50 ]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(413);
   plot(time, Bobs(:,3),'k',timesim, Bsim(:,3),'r','LineWidth',LW1); %Bz
   h.XTickLabel = [];
   h.Position(4) = Height;   
   ylabel('Bz [nT]');
   xlim([min(time) max(time)]);
   ylim([min(Bobs(:,3))-50 max(Bobs(:,3))+50 ]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(414);   
   plot(time,BobsStrength,'k',timesim,BsimStrength,'r','LineWidth',LW1); %|B|
   h.Position(4) = Height;   
   ylabel('B [nT]');
   xlim([min(time) max(time)]);
   ylim([min(BobsStrength)-50 max(BobsStrength)+50 ]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
end
