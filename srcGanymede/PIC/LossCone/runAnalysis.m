% Loss cone particle analysis
%
% 1. Read particle info from PIC
% 2. Read field info from PIC
% 3. Select particles in a given region; get B field at particle positions
% 4. Calculate pitch angles
% 5. Read field and topology info from MHD
% 6. Get Bmax = Bmax(theta,phi) on the surface
% 7. Find the corresponding field line footprint for each particle
% 8. Compute the surface B field strength for each particle
% 9. Get mirror ratio and loss cone angle for each particle
% 10.Select particles inside the loss cone.
% 11.Mapping particles onto the surface.
%
% 12.Calculate flux from bulk velocity and pressure.
%
%
% Hongyang Zhou, hyzhou@umich.edu 10/10/2018

clear; clc; close all
%%

% Get pitch angles for all particles
[angle,Bx_P,By_P,Bz_P,B_P,particle,weight] = get_pitch_angle;

% Get particles inside the loss cone
% [particle,vPerp,vPar,phi1,theta1] = ...
%    get_losscone(particle,angle,Bx_P,By_P,Bz_P,B_P,weight);

[particle,angle,Bsurf,Bx_P,By_P,Bz_P,B_P,theta1,phi1] = ...
   get_losscone(particle,angle,Bx_P,By_P,Bz_P,B_P,weight);


%%
calc_flux(particle,weight,Bx_P,By_P,Bz_P,B_P,Bsurf);

% Mapping particles onto the surface
%mapping_particle(particle,angle,phi1,theta1,Bsurf,B_P);

% Get bulk velocity and pressure
%p = get_pressure(particle,weight);

%%
figure
axesm('ortho','origin',[45 45]); 
axis off;
gridm on; 
framem on;
mlabel('equator')
plabel(0); 
plabel('fontweight','bold')

bins = 30;
phiMin = min(phi1); phiMax = max(phi1);
thetaMin = min(theta1); thetaMax = max(theta1);
dPhi = (phiMax - phiMin)/(bins - 1); 
dTheta = (thetaMax - thetaMin)/(bins - 1);

theta = thetaMin:dTheta:thetaMax;
phi = phiMin:dPhi:phiMax;
[Phi,Theta] = meshgrid(phi,theta);

[n,c] = hist3([phi1,theta1],'CdataMode','auto','Nbins',[30 30]);
h1 = surfm(Theta,Phi,n');
setm(gca,'Origin',[0 180 0])
plabel(0)
