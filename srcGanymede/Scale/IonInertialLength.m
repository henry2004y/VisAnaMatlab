% Estimation of plasma parameters near Ganymede.
%
% Ion inertial length near the magnetopause of Ganymede
%
% The smallest is about 0.06 RG, which is about 2 cells.
% 
% Hongyang Zhou, hyzhou@umich.edu 01/30/2018

clear; clc
%% Parameters
mu0      = 4*pi*1e-7;    %[H/m]
me       = 9.1094e-31;   %[kg] 
mp       = 1.6726*1e-27; %[kg]
epsilon0 = 8.8542*1e-12; %[F/m]
e        = 1.6022*1e-19; %[C]
c        = 2.9979*1e8;   %[m/s]
RG       = 2634000;      %[m]
Z        = 1;            % Charge state
mu       = 14;           % Average mass ratio

ipict = 1;

IonInertialsi = c / (sqrt(e^2/(epsilon0*mu*mp)*Z^2*4)*1e3*RG);
% 2.28*1e2/Z*sqrt(mu/56)*1e3/RG

%% Obtain the typical scales from simulation
filenamePC = ['~/Documents/research/Ganymede/data/EnergeticFlux/'...
   '3d_fluid_region0_0_t00000557_n00010710.out'];

[filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);

data = data.file1;
x = data.x(:,:,:,1); % [RG]
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

ne_  = strcmpi('rhos0',filehead.wnames);
ni_  = strcmpi('rhos1',filehead.wnames);
bx_  = strcmpi('bx',filehead.wnames);
by_  = strcmpi('by',filehead.wnames);
bz_  = strcmpi('bz',filehead.wnames);
uex_ = strcmpi('uxs0',filehead.wnames);
uey_ = strcmpi('uys0',filehead.wnames);
uez_ = strcmpi('uzs0',filehead.wnames);
uix_ = strcmpi('uxs1',filehead.wnames);
uiy_ = strcmpi('uys1',filehead.wnames);
uiz_ = strcmpi('uzs1',filehead.wnames);


Rhoe = data.w(:,:,:,ne_);    % [amu/cc]
Rhoi = data.w(:,:,:,ni_)/mu; % [amu/cc]
Bx   = data.w(:,:,:,bx_);    % [nT]
By   = data.w(:,:,:,by_);
Bz   = data.w(:,:,:,bz_);
Uex  = data.w(:,:,:,uex_);    % [km/s]
Uey  = data.w(:,:,:,uey_);
Uez  = data.w(:,:,:,uez_);
Uix  = data.w(:,:,:,uix_);    % [km/s]
Uiy  = data.w(:,:,:,uiy_);
Uiz  = data.w(:,:,:,uiz_);


x    = permute(x,[2 1 3]);
y    = permute(y,[2 1 3]);
z    = permute(z,[2 1 3]);
Rhoe = permute(Rhoe,[2 1 3]);
Rhoi = permute(Rhoi,[2 1 3]);
Bx   = permute(Bx,[2 1 3]);
By   = permute(By,[2 1 3]);
Bz   = permute(Bz,[2 1 3]);
Uex  = permute(Uex,[2 1 3]);
Uey  = permute(Uey,[2 1 3]);
Uez  = permute(Uez,[2 1 3]);
Uix  = permute(Uix,[2 1 3]);
Uiy  = permute(Uiy,[2 1 3]);
Uiz  = permute(Uiz,[2 1 3]);

% Ion plasma frequency
FreqIon = sqrt(e^2/(epsilon0*mu*mp)*Z^2*Rhoi)*1e3; %[rad/s]
% Ion inertial length, [m]
IonInertial = c / FreqIon;

% Ion inertial length normalized to Ganymede`s radius, [RG]
IonInertial = IonInertial / RG;

%sliceomatic(LambdaI);

%% Gyroradius 

B = sqrt(Bx.^2 + By.^2 + Bz.^2);
b = [Bx(:)'; By(:)'; Bz(:)'] ./ B(:)';

% Electron gyroradius, [RG]
UePar = abs(dot([Uex(:)'; Uey(:)'; Uez(:)'],b));
UePerp = sqrt(Uex(:)'.^2 + Uey(:)'.^2 + Uez(:)'.^2 - UePar.^2);
UePerp = reshape(UePerp,size(x));
Re = me.*UePerp*1e3./(e*B*1e-9) ./ RG;

% Ion gyroradius, [RG]
UiPar = abs(dot([Uix(:)'; Uiy(:)'; Uiz(:)'],b));
UiPerp = sqrt(Uix(:)'.^2 + Uiy(:)'.^2 + Uiz(:)'.^2 - UiPar.^2);
UiPerp = reshape(UiPerp,size(x));
Ri = mp*mu.*UiPerp*1e3./(e*B*1e-9) ./ RG;


%% Visualization
% Choose your cut
cut = 'y'; PlaneIndex = 64;
cut1 = squeeze(x(PlaneIndex,:,:));
cut2 = squeeze(z(PlaneIndex,:,:));
Ri = squeeze(Ri(PlaneIndex,:,:));

figure;
contourf(cut1,cut2,Ri,50,'Linestyle','none');
colorbar; axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
set(gca,'FontSize',16,'LineWidth',1.1)



