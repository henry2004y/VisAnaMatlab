% Loss cone particle analysis
%
% This script is used for analyzing particle information from MHD-EPIC
% simulations.
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
% 11.Calculate energy flux along the field line at the original locations
% 12.Map the fluxes onto the surface.
%
% Do I need to rewrite those parts that consider weights?
%
% Hongyang Zhou, hyzhou@umich.edu 10/10/2018

clear; clc; %close all
%% Get pitch angles for all particles

[angle,Bx_P,By_P,Bz_P,B_P,particle,weight] = getParticleInfo;

%% Get particles inside the loss cone

[particle,angle,Bsurf,Bx_P,By_P,Bz_P,B_P,theta,phi] = ...
   getLossCone(particle,angle,Bx_P,By_P,Bz_P,B_P,weight);

clearvars weight

%%
% Method 1
%calcEnergyFlux(particle,angle,true);
calcEnergyFlux(particle,angle);
% Method 2 (potentially wrong)
%calcEnergyFlux2(particle,angle,phi,theta,Bsurf,B_P);

% Test of mapping particles onto the surface
%mappingParticle(particle,angle,phi,theta,Bsurf,B_P);

return

%% Plot together

% Ion
figure

load('northIon.mat')

axesm('ortho','origin',[45 45]); 
axis off;
gridm on; 
framem on;
mlabel('equator')
surfm(theta,phi,flux); c = colorbar;

hold on

load('southIon.mat')

setm(gca,'Origin',[0 180 0],'Fontsize',12)
plabel(120); 
plabel('fontweight','bold')
surfm(theta,phi,flux); c = colorbar;
ylabel(c,'[W/m^2]'); c.FontSize = 14;
title('Ion Energetic Flux Density')

% Electron
figure

load('northElectron.mat')

axesm('ortho','origin',[45 45]); 
axis off;
gridm on; 
framem on;
mlabel('equator')
surfm(theta,phi,flux); c = colorbar;

hold on

load('southElectron.mat')

setm(gca,'Origin',[0 180 0],'Fontsize',12)
plabel(120); 
plabel('fontweight','bold')
surfm(theta,phi,flux); c = colorbar;
ylabel(c,'[W/m^2]'); c.FontSize = 14;
title('Electron Energetic Flux Density')

%%
% filePeak = '~/Documents/research/Ganymede/data/Aurora_peak_1998.csv';
% % PartI is mean, PartII is upper limit, PartIII is lower limit
% M = csvread(filePeak);
% 
% % Shift in longitude
% M(:,1) = M(:,1) - 90;
% 
% % Error bar
% %Error = M(15:28,2) - M(29:end,2);
% 
% %open('ElectronFlux.fig')
% open('IonFlux.fig')
% set(gcf, 'Position',  [100, 100, 650, 650])
% pos = plotm(M(1:14,2),M(1:14,1),'+r');
% 
% for i=1:14
%    linem(M(i+14:14:i+28,2),M(i+14:14:i+28,1),'k-','LineWidth',1)
% end
% 
% setm(gca,'FontSize',16)
% legend(pos,'Obs. Peak Emission')

%%
fileOval = '~/Documents/research/Ganymede/data/Oval Location Data.xlsx';

sheet = 1;

% 1998 north hemisphere
latlonNorth = xlsread(fileOval,sheet,'B6:C18');
latlonSouth = xlsread(fileOval,sheet,'B20:C32');

% plotm(51, 190,'*k')
% plotm(46,260,'*k')

% Inverse and shift the west longitude
latlonNorth(:,2) = 360 - (latlonNorth(:,2) - 90);
latlonSouth(:,2) = 360 - (latlonSouth(:,2) - 90);

%open('ElectronFlux.fig')
figure(4)
set(gcf, 'Position',  [100, 100, 650, 650])
pos = plotm(latlonNorth(:,1),latlonNorth(:,2),'+r');
pos = plotm(latlonSouth(:,1),latlonSouth(:,2),'+r');