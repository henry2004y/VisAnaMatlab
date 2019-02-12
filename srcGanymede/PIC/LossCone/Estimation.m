% Estimation

me = 9.10938356e-31; %[kg]
ve = 200e3; %[km/s]
q  = 6.242e+18; %[C]
N  = 1e3; %[cc]
we = N*1e6 * 0.5*me*ve^2 * q; %[eV]

%%
me = 9.10938356e-31; %[kg]

ne_surf = 100e6; %[#/m^3]

ve = 500e3; %[km/s]
% Planck constant
h = 6.626e-34; %[J s]
% wavelength of OI
lambda_OI = 1356e-10; %[m]

RG = 2634e3; %[m]
d =  0.3*RG; % assumed atomosphere column depth?

