% Estimation of plasma paramters near Ganymede.
%
% Ion inertial length near the magnetopause of Ganymede
%
% The smallest is about 0.06, which is about 2 cells.
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

LambdaIEsti = c / (sqrt(e^2/(epsilon0*mu*mp)*Z^2*4)*1e3*RG)
% 2.28*1e2/Z*sqrt(mu/56)*1e3/RG

%% Obtain the typical scales from simulation
filenamePC = '~/Ganymede/newPIC/run_G8_newPIC/3d_fluid_35.outs';

[filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);

data = data.file1;
x = data.x(:,:,:,1); % [RG]
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

Rhoe = data.w(:,:,:,1);    % [amu/cc]
Rhoi = data.w(:,:,:,2)/mu;    % [/cc]
Bx = data.w(:,:,:,3);      % [nT]
By = data.w(:,:,:,4);
Bz = data.w(:,:,:,5);
Uex = data.w(:,:,:,9);     % [km/s]
Uey = data.w(:,:,:,10);
Uez = data.w(:,:,:,11);


Rhoe = permute(Rhoe,[2 1 3]);
Rhoi = permute(Rhoi,[2 1 3]);
Bx = permute(Bx,[2 1 3]);
By = permute(By,[2 1 3]);
Bz = permute(Bz,[2 1 3]);
Uex = permute(Uex,[2 1 3]);
Uey = permute(Uey,[2 1 3]);
Uez = permute(Uez,[2 1 3]);

% Ion plasma frequency
OmegaI = sqrt(e^2/(epsilon0*mu*mp)*Z^2*Rhoi)*1e3; %[rad/s]
% Ion inertial length
LambdaI = c / OmegaI;

% Ion inertial length normalized to Ganymede`s radius
LambdaI = LambdaI / RG;

%sliceomatic(LambdaI);

% Electron gyroradius
B = sqrt(Bx.^2 + By.^2 + Bz.^2);
Ue = sqrt(Uex.^2 + Uey.^2 + Uez.^2);
Re = me.*Ue*1e3./(e*B*1e-9);

