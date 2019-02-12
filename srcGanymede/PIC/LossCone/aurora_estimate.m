% Estimate the contribution to Ganymede aurora emission.
%
% Brightness dependencies:
% 1. the number density of electrons (ne) with energies greater than 
% 14.3 eV penetrating in a column of neutral atmosphere at Ganymede;
% 2. the atmospheric column density of molecular oxygen at Ganymede, N(O2);
% 3. the collisional excitation rate of oxygen molecules by electrons 
% C(Te), and thus the electron temperatures (Te).
%
% Hongyang Zhou, hyzhou@umich.edu 01/05/2018

%% Parameters
kB = 1.38064852e-23; % [m^2*kg*s^{-2}*K^{-1}]
q  = 1.602e-19;      % [C]
me = 9.10938356e-31; % [kg]
mp = 1.6726219e-27;  % [kg]

a0 = 5.29e-9; % [cm], Bohr's radius

%% From Payan's paper, equation 3

% Energy to velocity conversion, the minimum speed required to excite O
Ethres = 14.3; % [eV]
vMin = sqrt(Ethres*q*2/me)/1e3; % [km/s]

% Collisional excitation rate of oxygen
C = 8.0104e-8*sigma*kB*Te/(pi*a0^2*sqrt(kB*Te))*exp(-Ethres/kB*Te);

%%
T = readtable('crosssection.csv','Delimiter',',',...
   'TreatAsEmpty',{'.','NA',''},'Format','%f%f%f%f')
