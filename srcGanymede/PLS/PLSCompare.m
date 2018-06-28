% Ganymede flyby simulation comparison along Galileo trajectory with PLS.
%
% Primarily used for G1 and G2 flybys.
%
% Hongyang Zhou, hyzhou@umich.edu 02/24/2018

clear;clc
%% Parameters
flyby = 2;   % [1,2,7,8,28,29]
DoPlot = 1;  % Plot output

%% Read observation data
flybyfile = strcat('Galileo_G',int2str(flyby),'_flyby_MAG.dat');
f = fullfile('~/Ganymede/GalileoData/galileomagdata',flybyfile);
[~,data] = read_log_data(f);

time = datetime(data(:,1:6));
xyz  = data(:,7:9);
Bobs = data(:,10:12);
BobsMag = sqrt(Bobs(:,1).^2+Bobs(:,2).^2+Bobs(:,3).^2);
%filename = '~/Ganymede/newPIC/G2/box_traj_steady.out';
%filename = '~/Ganymede/newPIC/G2/box_var_4_t00000800_n00221232.out';
filename='~/Ganymede/newPIC/G2/box_fixed_test.out';


npict = 1; % Remember to change this for different runs!
[filehead,data] = read_data(filename,'npict',1);

% Interpolate simulation data to observation data
data = data.file1;
nx = filehead.nx; nw = filehead.nw;

x = data.x(:,:,:,1);
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

w = data.w;
% Find the correct index of variables
bx_ = strcmpi('bx',filehead.wnames);
by_ = strcmpi('by',filehead.wnames);
bz_ = strcmpi('bz',filehead.wnames);
ux_ = strcmpi('ux',filehead.wnames);
uy_ = strcmpi('uy',filehead.wnames);
uz_ = strcmpi('uz',filehead.wnames);

bx = data.w(:,:,:,bx_);
by = data.w(:,:,:,by_);
bz = data.w(:,:,:,bz_);
ux = data.w(:,:,:,ux_);
uy = data.w(:,:,:,uy_);
uz = data.w(:,:,:,uz_);

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

% gridded interpolation, much faster than scattered interpolation
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
   
%% Compare with PLS data
f = '~/Ganymede/newPIC/G2/G2_Frank1997.dat';
[head,data] = read_log_data(f);

timePLS = datetime(data(:,1:6));
Uobs = data(:,7:9);
UobsMag = sqrt(Uobs(:,1).^2+Uobs(:,2).^2+Uobs(:,3).^2);

% Assuming heavy ions mi=14 inside the magnetosphere
mi = 14;
UObsIon = Uobs(18:43,:) / sqrt(mi);
UObsIonMag = sqrt(UObsIon(:,1).^2+UObsIon(:,2).^2+UObsIon(:,3).^2);

%% Visualization
if DoPlot
   figure('Position', [100, 100, 1100, 800]);
   h = subplot(411); LW1 = 2.5; LW2 = 1.5; FS=16; Height = 0.18;
   plot(time,Usim(:,1),'r','LineWidth',LW1); %Ux 
   hold on
   plot(timePLS,Uobs(:,1),'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
   plot(timePLS(18:43),UObsIon(:,1),'.-k','MarkerSize',16)
   hline(0,'--k')
   h.XTickLabel = []; 
   h.Position(4) = Height;
   ylabel('Ux [km/s]');
   xlim([min(time) max(time)]);
   %ylim([min(Bobs(:,1))-50 max(Bobs(:,1))+50 ]);
   %title({'(a)',strcat('Galileo G',int2str(flyby),' Flyby Magnetic field')})
   title(strcat('Galileo G',int2str(flyby),' Flyby Magnetic field'))
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(412);
   plot(time,Usim(:,2),'r','LineWidth',LW1); %Uy
   hold on
   plot(timePLS,Uobs(:,2),'.-','MarkerSize',16,'Color',[0.5 0.5 0.5])
   plot(timePLS(18:43),UObsIon(:,2),'.-k','MarkerSize',16)
   hline(0,'--k')
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
   hline(0,'--k')
   
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

%% Parallel and perpendicular components of the bulk flow 
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
hline(0,'--k')

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



