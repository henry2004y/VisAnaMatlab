clear; clc
% Inner boundary setup

k = 1.38e-23; % [J/k]
g = 1.428;    % [m/s^2]
mp= 1.6726e-27; % [kg]
e = 1.6022e-19; % [C]

%% Upstream condition

% G2 flyby
piG2 = (1e6*k*1e6+1e6*k*5e6)*1e9

% upstream thermal pressure
piUp = 3.6; %[nPa]
peUp = 0.2; %[nPa]
nUp = 4;   %[cm^-3]


TiUp = piUp*1e-9 / (nUp*1e6 * k)
TeUp = peUp*1e-9 / (nUp*1e6 * k)
% upstream magnetic pressure
PbUp = sqrt(10^2+6^2+86^2)*(1e-9)^2/(8*pi*1e-7)*1e9


%% Xianzhe`s parameter
% n=39cm^-3, kT=20eV, average ion mass = 14
Told   = 20*e/k  %[K]
HscaleOld = 20*e/(14*mp*g)*1e-3 %[km]
pold   = 550*1e6/14*k*Told*1e9 %[nPa]

%% Test parameter

% upstream
TiUp = piUp*1e-9 / (nUp*1e6 * k)

n = 100; % [cm^-3]
T = 4e4;   % [K]

Hscale = T*k/(14*mp*g)*1e-3 %[km]
p = n*1e6*k*T*1e9

%%
T = 800; %[K]
n = 39; %[cc], rho=550
Hscale = T*k/(14*mp*g)*1e-3 %[km]
p = n*1e6*k*T*1e9

%%
T = 800; %[K]
n = 100; %[cc], rho=1400
Hscale = T*k/(14*mp*g)*1e-3 %[km]
p = n*1e6*k*T*1e9

%%
T = 800; %[K]
n = 100/14; %[cc], rho=100
Hscale = T*k/(14*mp*g)*1e-3 %[km]
p = n*1e6*k*T*1e9

%%
Rg = 2634000;
mu0 = 4*pi*1e-7;
eta = 6e11 * mu0;

DiffusionTime = (0.5*Rg)^2 / eta;

