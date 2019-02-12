% Estimation of surface electron density from PIC simulation.

clear; clc; close all
%% Get all particle info

[angle,Bx_P,By_P,Bz_P,B_P,particle,weight] = getParticleInfo;

vx_ = Parameters.vx_;
vy_ = Parameters.vy_;
vz_ = Parameters.vz_;
box = Parameters.Region;

q         = Parameters.q;
kB        = Parameters.kB;
me        = Parameters.me;
mp        = Parameters.mp;
mi        = Parameters.mi;
No2SiMass = Parameters.No2SiMass;
Rg        = Parameters.Rg;

uBulk = mean(particle(vx_:vz_,:),2);

Volume = (box(2)-box(1))*(box(4)-box(3))*(box(6)-box(5))*Rg^3; %[m^3]

rhoMean = No2SiMass * sum(abs(weight)) / Volume / mp / 1e6;

uxMean = mean(particle(vx_,:));
uyMean = mean(particle(vy_,:));
uzMean = mean(particle(vz_,:));

% This temperature estimation is a little bit high.
TMean = 0.5*sum(sum(particle(4:6,:).^2,1).*abs(weight'))*No2SiMass / ...
   (sum(abs(weight')./0.14)*No2SiMass/mp) / q; % [eV]

PMean = 0.5*sum(sum(particle(4:6,:).^2,1).*abs(weight'))*No2SiMass / ...
   Volume *1e9; % [nPa]

%%
[particle,angle,Bsurf,Bx_P,By_P,Bz_P,B_P,theta,phi] = ...
   getLossCone(particle,angle,Bx_P,By_P,Bz_P,B_P,weight);

clearvars weight

%% Obtain topology for the mesh
Dir = Parameters.Dir;
fnameGM = Parameters.fnameGM;

[filehead,data] = read_data(fullfile(Dir,fnameGM),'verbose',false);
data = data.file1;

if strcmp(Parameters.Hemisphere,'north')
   theta_ = strcmpi('theta1',filehead.wnames);
   phi_ = strcmpi('phi1',filehead.wnames);
else
   theta_ = strcmpi('theta2',filehead.wnames);
   phi_ = strcmpi('phi2',filehead.wnames);
end

xGM = data.x(:,:,:,1);       % [Rg]
yGM = data.x(:,:,:,2);       % [Rg]
zGM = data.x(:,:,:,3);       % [Rg]

theta = data.w(:,:,:,theta_);
phi = data.w(:,:,:,phi_);

Ftheta1 = griddedInterpolant(xGM,yGM,zGM,theta);
Fphi1 = griddedInterpolant(xGM,yGM,zGM,phi);

clearvars Dir fnameGM theta1_ phi1_ data xGM yGM zGM theta1 phi1

xMin = min(particle(1,:)); xMax = max(particle(1,:));
yMin = min(particle(2,:)); yMax = max(particle(2,:));
zMin = min(particle(3,:)); zMax = max(particle(3,:));

bins = Parameters.bins;

y = linspace(yMin,yMax,bins);
z = linspace(zMin,zMax,bins);
[Y,Z] = ndgrid(y,z);
dy = (yMax - yMin)/(bins - 1); 
dz = (zMax - zMin)/(bins - 1);
subs = [round((particle(2,:)'-yMin)/dy)+1 ...
        round((particle(3,:)'-zMin)/dz)+1];

% Volume, [m^3]
Volume = (xMax - xMin)*dy*dz * Parameters.Rg^3;

rho = accumarray(subs, abs(particle(7,:))).* No2SiMass / mp /1e6 ./ Volume;

ne = rho ./ 0.14;

%% Mapping to the surface
X = mean(particle(1,:)) * ones(size(Y));
theta = Ftheta1(X,Y,Z);
phi   = Fphi1(X,Y,Z);

figure
surf(phi,theta,ne); colorbar

% Get B field at the orginal locations
[xF,yF,zF,Bx,By,Bz] = getField;
Fx = griddedInterpolant(xF,yF,zF,Bx);
Fy = griddedInterpolant(xF,yF,zF,By);
Fz = griddedInterpolant(xF,yF,zF,Bz);

Bx = Fx(X,Y,Z); By = Fy(X,Y,Z); Bz = Fz(X,Y,Z);
B  = sqrt(Bx.^2 + By.^2 + Bz.^2);

[FBxSurf,FBySurf,FBzSurf] = getBsurface(false);
BxSurf = FBxSurf(phi,theta);
BySurf = FBySurf(phi,theta);
BzSurf = FBzSurf(phi,theta);
Bsurf = sqrt(BxSurf.^2 + BySurf.^2 + BzSurf.^2);

% This may not be valid! n/B ~= const.
neSurf = ne .* Bsurf ./ B;

% figure
% surf(phi,theta,Bsurf./B); colorbar
% 
figure
surf(phi,theta,neSurf); colorbar
title('Ne surf')

%% Brightness estimation
NO2 = 4e14; % [cm^{-2}]
C = 5e-8;
Brightness = neSurf*C*NO2*1e-6; % [R]

figure
surf(phi,theta,Brightness); colorbar
title('Brightness [R]')


% Secondary excitation
% As being discussed in Eviator's paper, collisions are rare in Ganymede's
% atmosphere because of the low column density. So the secondary electrons
% are not important here!
C2 = 1e-8;
Brightness2 = neSurf*TMean/35 * C2 * NO2 *1e-6;

figure
surf(phi,theta,Brightness2); colorbar
title('Brightness2 [R]')