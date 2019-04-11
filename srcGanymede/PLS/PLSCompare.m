% Ganymede flyby simulation comparison along Galileo trajectory.
% Galileo PLS data are only available for G1 and G2 flyby.
%
% Data sources:
% 1. Paper "New Results From Galileo's First Flyby of Ganymede:
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
% I don't have a MHD-EPIC simulation for G1, G2 or G7, where PLS or EPD
% data are presented in literature. The Hall MHD results show that in
% general it has stronger reconnection signature than ideal MHD even on the
% flank of the magnetopause, which may not be observed by some
% postprocessed plasma data. 
% It is useful to: 
% 1. perform a PIC simulation using a PIC box on the side;
% 2. compare ideal/Hall/PIC simulation.
%
% Hongyang Zhou, hyzhou@umich.edu 05/13/2018

clear;clc;%close all
%% Parameters
kb           = 1.38e-23;
IonMassRatio = 14;
PLSData      = 2;
flyby        = 2;   % [1,2]
PLS          = 'S5'; % S1 - S8
DoPlot       = false;  % Plot output
DoPlotObs    = false;
DoPlotParPerp= true;
npict        = 1; % Remember to change this for different simulation runs!

if PLSData == 1
   IndexIn_ = 15:46;
else
   IndexIn_ = 18:43;
end

% Figure parameters
LW1 = 2.5; LW2 = 1.5; FS=16; Height = 0.18;
   
%% Read Galileo trajectory data
flybyfile = strcat('Galileo_G',int2str(flyby),'_flyby_MAG.dat');
fileMag = fullfile('/Users/hyzhou/Ganymede/GalileoData/galileomagdata',...
   flybyfile);
[~,data] = read_log_data(fileMag);

time = datetime(data(:,1:6));
xyz  = data(:,7:9);
Bobs = data(:,10:12);
BobsMag = sqrt(Bobs(:,1).^2+Bobs(:,2).^2+Bobs(:,3).^2);

%% Read simulation data
filename = '~/Ganymede/newPIC/G2/box_traj_steady.out';
%filename = '~/Ganymede/MOP2018/runG7_steady/GM/box*.outs'; % Maybe not fully covered?
%filename='~/SWMF/SWMF/GM/BATSRUS/run_test/RESULTS/run_G1_test3/GM/box*';
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

%% Read PLS data

if PLSData == 1
% Import the data
[~, ~, raw] = xlsread(...
   '/Users/hyzhou/Ganymede/GalileoData/GalileoPLS_G1G2.xlsx',PLS);
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
Obs.V  = sqrt(Obs.Vx.^2 + Obs.Vy.^2 + Obs.Vz.^2);

% Clear temporary variables
clearvars data raw stringVectors

% Plot Obs Data
if DoPlotObs
   figure('Position', [100, 100, 1100, 800]);
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
   plot(time,Usim(:,1),'r','LineWidth',LW1); %Ux
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
   plot(time,Usim(:,2),'r','LineWidth',LW1); %Uy
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
   plot(time, Usim(:,3),'r','LineWidth',LW1); %Uz
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
   plot(time,UsimMag,'r','LineWidth',LW1); %|U|
   hold on
   plot(Obs.time,Obs.V,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5]);
   %h.Position(4) = Height;
   ylabel('U [km/s]');
   xlim([min(time) max(time)]);
   %ylim([0 150]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
end


elseif PLSData == 2
   f = '~/Ganymede/newPIC/G2/G2_Frank1997.dat';
   [head,data] = read_log_data(f);
   
   Obs.time = datetime(data(:,1:6));
   Obs.Vx   = data(:,7);
   Obs.Vy   = data(:,8);
   Obs.Vz   = data(:,9);
   Obs.V    = sqrt(Obs.Vx.^2+Obs.Vy.^2+Obs.Vz.^2);
   
   % Assuming heavy ions mi=14 inside the magnetosphere
   Obs.Vix = Obs.Vx(IndexIn_,:) / sqrt(IonMassRatio);
   Obs.Viy = Obs.Vy(IndexIn_,:) / sqrt(IonMassRatio);
   Obs.Viz = Obs.Vz(IndexIn_,:) / sqrt(IonMassRatio);
   Obs.Vi = sqrt(Obs.Vix.^2+Obs.Viy.^2+Obs.Viz.^2);
   
   if DoPlotObs
      figure('Position', [100, 100, 1100, 800]);
      h = subplot(411); LW1 = 2.5; LW2 = 1.5; FS=16; Height = 0.18;
      plot(time,Usim(:,1),'r','LineWidth',LW1); %Ux
      hold on
      plot(Obs.time,Obs.Vx,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
      plot(Obs.time(IndexIn_),Obs.Vix,'.-k','MarkerSize',16)
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
      plot(Obs.time,Obs.Vy,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
      plot(Obs.time(IndexIn_),Obs.Viy,'.-k','MarkerSize',16)
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
      plot(Obs.time,Obs.Vz,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
      plot(Obs.time(IndexIn_),Obs.Viz,'.-k','MarkerSize',16)
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
      plot(Obs.time,Obs.V,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5]);
      plot(Obs.time(IndexIn_),Obs.Vi,'.-k','MarkerSize',16)
      h.Position(4) = Height;
      ylabel('U [km/s]');
      xlim([min(time) max(time)]);
      %ylim([min(BobsStrength)-50 max(BobsStrength)+50 ]);
      set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   end
   
else
   error('unknown PLS data source!') 
end


%% Visualization
if DoPlot
   figure('Position', [100, 100, 1100, 800]);
   
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
   plot(Obs.time,Obs.V,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5]);
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


BxObs = interp1(datenum(time),Bobs(:,1),datenum(Obs.time));
ByObs = interp1(datenum(time),Bobs(:,2),datenum(Obs.time));
BzObs = interp1(datenum(time),Bobs(:,3),datenum(Obs.time));
BobsMag = sqrt(BxObs.^2 + ByObs.^2 + BzObs.^2);

UparObs = (Obs.Vx.*BxObs + Obs.Vy.*ByObs + Obs.Vz.*BzObs)./BobsMag;
UperpObs = sqrt(Obs.V.^2 - UparObs.^2);

UparObsIon = UparObs ./ sqrt(IonMassRatio);
UperpObsIon = UperpObs ./ sqrt(IonMassRatio);

figure
Height = 0.25;
h = subplot(311);
plot(time,Upar,'r','LineWidth',LW1);
hold on
plot(Obs.time,UparObs,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5]);
plot(Obs.time,UparObsIon,'.-k','MarkerSize',16)
yline(0,'--k')

h.XTickLabel = [];
ylabel('$U_\parallel$ [km/s]','Interpreter','Latex')
h.Position(4) = Height;
set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);

h = subplot(312);
plot(time,Uperp,'r','LineWidth',LW1);
hold on
plot(Obs.time,UperpObs,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5]);
plot(Obs.time,UperpObsIon,'.-k','MarkerSize',16)
h.XTickLabel = [];
ylabel('U_\perp [km/s]')
h.Position(4) = Height;
set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);

h = subplot(313);
plot(time,UsimMag,'r','LineWidth',LW1); %|U|
hold on
plot(Obs.time,Obs.V,'.-','MarkerSize',16,'Color',[0.5 0.5 0.5]);
plot(Obs.time,Obs.V./sqrt(IonMassRatio),'.-k','MarkerSize',16)
h.Position(4) = Height;
ylabel('U [km/s]');
xlim([min(time) max(time)]);

set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);

