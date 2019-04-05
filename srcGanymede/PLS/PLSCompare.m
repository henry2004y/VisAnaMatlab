% Ganymede flyby simulation comparison along Galileo trajectory.
% Galileo PLS data are only available for G1 and G2 flyby.
%
% Data sources:
% 1. Paper "New Results From Galileo?s First Flyby of Ganymede:
% Reconnection-Driven Flows at the Low-Latitude Magnetopause Boundary,
% Crossing the Cusp, and Icy Ionospheric Escape"
%
% PLS data set:
% S0. Single Spin G01 moments assuming mass per charge (Mi/Qi) = 16
% S1. G01 Flyby Assuming mass per charge (Mi/Qi) = 1
% S2. G01 Flyby Assuming mass per charge (Mi/Qi) = 2
% S3. G01 Flyby Assuming mass per charge (Mi/Qi) = 16
% S4. G01 Flyby Assuming mass per charge (Mi/Qi) = 32
% S5. G02 Flyby Assuming mass per charge (Mi/Qi) = 1
% S6. G02 Flyby Assuming mass per charge (Mi/Qi) = 2
% S7. G02 Flyby Assuming mass per charge (Mi/Qi) = 16
% S8. G02 Flyby Assuming mass per charge (Mi/Qi) = 32
%
% 2. Internal communication, obtained from L. A. Frank, G2 flyby only 
%
% Hongyang Zhou, hyzhou@umich.edu 05/13/2018

clear;clc
%% Parameters
kb           = 1.38e-23;
IonMassRatio = 14;
PLSData      = 1;
flyby        = 2;   % [1,2]
PLS          = 'S7'; % S1 - S8
DoPlot       = false;  % Plot output
DoPlotObs    = 1;
DoPlotParPerp= false;
npict        = 1; % Remember to change this for different simulation runs!

%% Read Galileo trajectory data
flybyfile = strcat('Galileo_G',int2str(flyby),'_flyby_MAG.dat');
fileMag = fullfile('../Galileo',flybyfile);
[~,data] = read_log_data(fileMag);

time = datetime(data(:,1:6));
xyz  = data(:,7:9);
Bobs = data(:,10:12);
BobsMag = sqrt(Bobs(:,1).^2+Bobs(:,2).^2+Bobs(:,3).^2);

%% Read PLS data

if PLSData == 1
% Import the data
[~, ~, raw] = xlsread(...
   '/Users/hyzhou/Documents/research/Ganymede/data/GalileoPLS_G1G2',PLS);
raw = raw(3:end,:);
stringVectors = string(raw(:,1));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,[2,3,4,5,6,7,8,9]);

% Create output variable
data = reshape([raw{:}],size(raw));

% Create table
Obs = table;

% Allocate imported array to column variable names
Obs.time = datetime(stringVectors(:,1),...
   'InputFormat','uuuu-DDD''//''HH:mm:ss');
Obs.X = data(:,1);
Obs.Y = data(:,2);
Obs.Z = data(:,3);
Obs.N = data(:,4);
Obs.T = data(:,5)/1.16045221e4; % [K] --> [eV]
Obs.Vx = data(:,6);
Obs.Vy = data(:,7);
Obs.Vz = data(:,8);

% Clear temporary variables
clearvars data raw stringVectors


% Plot Obs Data
if DoPlotObs
   figure('Position', [100, 100, 1100, 800]);
   LW1 = 2.5; LW2 = 1.5; FS=16; %Height = 0.14;
   title(strcat('Galileo G',int2str(flyby),' Flyby'))
   h = subplot(611); 

   hold on
   plot(Obs.time,Obs.N,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
   set(gca, 'YScale', 'log');
   h.XTickLabel = [];
   %h.Position(4) = Height;
   ylabel('Density');
   xlim([min(time) max(time)]);
   %ylim([1 3e2]);
   %title({'(a)',strcat('Galileo G',int2str(flyby),' Flyby Magnetic field')})
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(612); 

   hold on
   plot(Obs.time,Obs.T,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
   set(gca, 'YScale', 'log');
   h.XTickLabel = [];
   %h.Position(4) = Height;
   ylabel('Temperature');
   xlim([min(time) max(time)]);
   %ylim([1 1e4]);
   %title({'(a)',strcat('Galileo G',int2str(flyby),' Flyby Magnetic field')})
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(613); 
   %plot(time,Usim(:,1),'r','LineWidth',LW1); %Ux
   hold on
   plot(Obs.time,Obs.Vx,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
   yline(0,'--k')
   h.XTickLabel = [];
   %h.Position(4) = Height;
   ylabel('Ux [km/s]');
   xlim([min(time) max(time)]);
   %ylim([-50 150]);
   %title({'(a)',strcat('Galileo G',int2str(flyby),' Flyby Magnetic field')})
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(614);
   %plot(time,Usim(:,2),'r','LineWidth',LW1); %Uy
   hold on
   plot(Obs.time,Obs.Vy,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
   yline(0,'--k')
   h.XTickLabel = [];
   %h.Position(4) = Height;
   ylabel('Uy [km/s]');
   xlim([min(time) max(time)]);
   %ylim([-50 50]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(615);
   %plot(time, Usim(:,3),'r','LineWidth',LW1); %Uz
   hold on
   plot(Obs.time,Obs.Vz,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
   yline(0,'--k')
   
   h.XTickLabel = [];
   %h.Position(4) = Height;
   ylabel('Uz [km/s]');
   xlim([min(time) max(time)]);
   %ylim([-50 50]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(616);
   %plot(time,UsimMag,'r','LineWidth',LW1); %|U|
   hold on
   plot(Obs.time,sqrt(Obs.Vx.^2+Obs.Vy.^2+Obs.Vz.^2),...
      '.-','MarkerSize',16,'Color',[0.5 0.5 0.5]);
   %h.Position(4) = Height;
   ylabel('U [km/s]');
   xlim([min(time) max(time)]);
   %ylim([0 150]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
end


elseif PLSData == 2
   f = '~/Ganymede/newPIC/G2/G2_Frank1997.dat';
   [head,data] = read_log_data(f);
   
   timePLS = datetime(data(:,1:6));
   Uobs = data(:,7:9);
   UobsMag = sqrt(Uobs(:,1).^2+Uobs(:,2).^2+Uobs(:,3).^2);
   
   % Assuming heavy ions mi=14 inside the magnetosphere
   UObsIon = Uobs(18:43,:) / sqrt(IonMassRatio);
   UObsIonMag = sqrt(UObsIon(:,1).^2+UObsIon(:,2).^2+UObsIon(:,3).^2);
   
   if DoPlotObs
      figure('Position', [100, 100, 1100, 800]);
      h = subplot(411); LW1 = 2.5; LW2 = 1.5; FS=16; Height = 0.18;
      plot(time,Usim(:,1),'r','LineWidth',LW1); %Ux
      hold on
      plot(timePLS,Uobs(:,1),'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
      plot(timePLS(18:43),UObsIon(:,1),'.-k','MarkerSize',16)
      yline(0,'--k')
      h.XTickLabel = [];
      h.Position(4) = Height;
      ylabel('Ux [km/s]');
      xlim([min(time) max(time)]);
      title(strcat('Galileo G',int2str(flyby),' Flyby'))
      set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
      h = subplot(412);
      plot(time,Usim(:,2),'r','LineWidth',LW1); %Uy
      hold on
      plot(timePLS,Uobs(:,2),'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
      plot(timePLS(18:43),UObsIon(:,2),'.-k','MarkerSize',16)
      yline(0,'--k')
      h.XTickLabel = [];
      h.Position(4) = Height;
      ylabel('Uy [km/s]');
      xlim([min(time) max(time)]);
      %ylim([min(Bobs(:,2))-50 max(Bobs(:,2))+50 ]);
      set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
      h = subplot(413);
      plot(time, Usim(:,3),'r','LineWidth',LW1); %Uz
      hold on
      plot(timePLS,Uobs(:,3),'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
      plot(timePLS(18:43),UObsIon(:,3),'.-k','MarkerSize',16)
      yline(0,'--k')
      
      h.XTickLabel = [];
      h.Position(4) = Height;
      ylabel('Uz [km/s]');
      xlim([min(time) max(time)]);
      %ylim([min(Bobs(:,3))-50 max(Bobs(:,3))+50 ]);
      set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
      h = subplot(414);
      plot(time,UsimMag,'r','LineWidth',LW1); %|U|
      hold on
      plot(timePLS,UobsMag(:),'.-','MarkerSize',16,'Color',[0.5 0.5 0.5]);
      plot(timePLS(18:43),UObsIonMag,'.-k','MarkerSize',16)
      h.Position(4) = Height;
      ylabel('U [km/s]');
      xlim([min(time) max(time)]);
      %ylim([min(BobsStrength)-50 max(BobsStrength)+50 ]);
      set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   end
   
else
   error('unknown PLS data source!') 
end

%% Read simulation data
filename = '~/Ganymede/newPIC/G2/box_traj_steady.out';
filename='~/SWMF/SWMF/GM/BATSRUS/run_test/RESULTS/run_G1_test3/GM/box*';
filename='~/Documents/research/Ganymede/SteadyRun/run_G1_test1/GM/box*';
[filehead,data] = read_data(filename,'npict',npict);

% Interpolate simulation data to observation data
data = data.file1;
%nx = filehead.nx; nw = filehead.nw;

x = data.x(:,:,:,1);
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

% Find the correct index of variables
bx_ = strcmpi('bx',filehead.wnames);
by_ = strcmpi('by',filehead.wnames);
bz_ = strcmpi('bz',filehead.wnames);
ux_ = strcmpi('ux',filehead.wnames);
uy_ = strcmpi('uy',filehead.wnames);
uz_ = strcmpi('uz',filehead.wnames);
rho_= strcmpi('rho',filehead.wnames);
p_  = strcmpi('p',filehead.wnames);

bx = data.w(:,:,:,bx_);
by = data.w(:,:,:,by_);
bz = data.w(:,:,:,bz_);
ux = data.w(:,:,:,ux_);
uy = data.w(:,:,:,uy_);
uz = data.w(:,:,:,uz_);
rho= data.w(:,:,:,rho_);
p  = data.w(:,:,:,p_);

% From ndgrid to meshgrid format
x  = permute(x,[2 1 3]);
y  = permute(y,[2 1 3]);
z  = permute(z,[2 1 3]);
bx = permute(bx,[2 1 3]);
by = permute(by,[2 1 3]);
bz = permute(bz,[2 1 3]);
ux = permute(ux,[2 1 3]);
uy = permute(uy,[2 1 3]);
uz = permute(uz,[2 1 3]);
rho= permute(rho,[2 1 3]);
p  = permute(p,[2 1 3]);

% Interpolate simulation data to observation data
% Gridded interpolation, much faster than scattered interpolation
Bsim = Inf(size(xyz,1),3);
Bsim(:,1) = interp3(x,y,z,bx,xyz(:,1),xyz(:,2),xyz(:,3));
Bsim(:,2) = interp3(x,y,z,by,xyz(:,1),xyz(:,2),xyz(:,3));
Bsim(:,3) = interp3(x,y,z,bz,xyz(:,1),xyz(:,2),xyz(:,3));
Usim = Inf(size(xyz,1),3);
Usim(:,1) = interp3(x,y,z,ux,xyz(:,1),xyz(:,2),xyz(:,3));
Usim(:,2) = interp3(x,y,z,uy,xyz(:,1),xyz(:,2),xyz(:,3));
Usim(:,3) = interp3(x,y,z,uz,xyz(:,1),xyz(:,2),xyz(:,3));

BsimMag = sqrt(Bsim(:,1).^2+Bsim(:,2).^2+Bsim(:,3).^2);
UsimMag = sqrt(Usim(:,1).^2+Usim(:,2).^2+Usim(:,3).^2);

rho = interp3(x,y,z,rho,xyz(:,1),xyz(:,2),xyz(:,3))/IonMassRatio;
p   = interp3(x,y,z,p,xyz(:,1),xyz(:,2),xyz(:,3));
T   = p*1e-9./(rho*1e6*kb)/1.16045221e4; %[K] --> [eV]

%% Visualization
if DoPlot
   figure('Position', [100, 100, 1100, 800]);
   LW1 = 2.5; LW2 = 1.5; FS=16; Height = 0.18;
   
   h = subplot(611);   
   plot(time,rho,'r','LineWidth',LW1); % rho
   title(strcat('Galileo G',int2str(flyby),' Flyby'))
   hold on
   plot(Obs.time,Obs.N,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
   set(gca, 'YScale', 'log');
   h.XTickLabel = [];
   %h.Position(4) = Height;
   ylabel('Density');
   xlim([min(time) max(time)]);
   %ylim([1 4e2]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   
   h = subplot(612); 
   plot(time,T,'r','LineWidth',LW1); % rho
   hold on
   plot(Obs.time,Obs.T,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
   set(gca, 'YScale', 'log');
   h.XTickLabel = [];
   %h.Position(4) = Height;
   ylabel('Temperature');
   xlim([min(time) max(time)]);
   %ylim([1 1e4]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);

   h = subplot(613); 
   plot(time,Usim(:,1),'r','LineWidth',LW1); %Ux 
   hold on
   plot(Obs.time,Obs.Vx,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
   yline(0,'--k')
   h.XTickLabel = []; 
   %h.Position(4) = Height;
   ylabel('Ux [km/s]');
   xlim([min(time) max(time)]);
   %ylim([min(Obs.Vx)-50 max(Obs.Vx)+50 ]);
   %ylim([-50 150]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   
   h = subplot(614);
   plot(time,Usim(:,2),'r','LineWidth',LW1); %Uy
   hold on
   plot(Obs.time,Obs.Vy,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
   yline(0,'--k')
   h.XTickLabel = [];
   %h.Position(4) = Height;   
   ylabel('Uy [km/s]');
   xlim([min(time) max(time)]);
   %ylim([min(Obs.Vy)-50 max(Obs.Vy)+50 ]);
   %ylim([-50 50]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   
   h = subplot(615);
   plot(time, Usim(:,3),'r','LineWidth',LW1); %Uz
   hold on
   plot(Obs.time,Obs.Vz,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
   yline(0,'--k')  
   h.XTickLabel = [];
   %h.Position(4) = Height;   
   ylabel('Uz [km/s]');
   xlim([min(time) max(time)]);
   %ylim([min(Obs.Vz)-50 max(Obs.Vz)+50 ]);
   %ylim([-50 50]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   
   h = subplot(616);   
   plot(time,UsimMag,'r','LineWidth',LW1); %|U|
   hold on
   plot(Obs.time,sqrt(Obs.Vx.^2+Obs.Vy.^2+Obs.Vz.^2),...
      '.-','MarkerSize',16,'Color',[0.5 0.5 0.5]);
   %h.Position(4) = Height;   
   ylabel('U [km/s]');
   xlim([min(time) max(time)]);
   %ylim([min(BobsStrength)-50 max(BobsStrength)+50 ]);
   %ylim([0 150]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
end

%% Parallel and perpendicular components of the bulk flow
if ~DoPlotParPerp, return, end

Upar = (Usim(:,1).*Bsim(:,1) + Usim(:,2).*Bsim(:,2) + ...
   Usim(:,3).*Bsim(:,3))./BsimMag;
Uperp = sqrt(UsimMag.^2 - Upar.^2);


BxObs = interp1(datenum(time),Bobs(:,1),datenum(timePLS));
ByObs = interp1(datenum(time),Bobs(:,2),datenum(timePLS));
BzObs = interp1(datenum(time),Bobs(:,3),datenum(timePLS));
BobsMag = sqrt(BxObs.^2 + ByObs.^2 + BzObs.^2);

UparObs = (Uobs(:,1).*BxObs + Uobs(:,2).*ByObs + ...
   Uobs(:,3).*BzObs)./BobsMag;
UperpObs = sqrt(UobsMag.^2 - UparObs.^2);

UparObsIon = UparObs(18:43)./ sqrt(mi);
UperpObsIon = UperpObs(18:43) ./ sqrt(mi);

figure
Height = 0.25;
h = subplot(311);
plot(time,Upar,'r','LineWidth',LW1);
hold on
plot(timePLS,UparObs,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5]);
plot(timePLS(18:43),UparObsIon,'.-k','MarkerSize',16)
yline(0,'--k')

h.XTickLabel = [];
ylabel('$U_\parallel$ [km/s]','Interpreter','Latex')
h.Position(4) = Height;
set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);

h = subplot(312);
plot(time,Uperp,'r','LineWidth',LW1);
hold on
plot(timePLS,UperpObs,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5]);
plot(timePLS(18:43),UperpObsIon,'.-k','MarkerSize',16)
h.XTickLabel = [];
ylabel('U_\perp [km/s]')
h.Position(4) = Height;
set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);

h = subplot(313);
plot(time,UsimMag,'r','LineWidth',LW1); %|U|
hold on
plot(timePLS,UobsMag(:),'.-','MarkerSize',16,'Color',[0.5 0.5 0.5]);
plot(timePLS(18:43),UObsIonMag,'.-k','MarkerSize',16)
h.Position(4) = Height;
ylabel('U [km/s]');
xlim([min(time) max(time)]);

set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);

