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
% 12.Mapping the fluxes onto the surface.
%
%
% Hongyang Zhou, hyzhou@umich.edu 10/10/2018

clear; clc; close all
%%

% Get pitch angles for all particles
[angle,Bx_P,By_P,Bz_P,B_P,particle,weight] = getParticleInfo;

% Get particles inside the loss cone
[particle,angle,Bsurf,Bx_P,By_P,Bz_P,B_P,theta1,phi1] = ...
   getLossCone(particle,angle,Bx_P,By_P,Bz_P,B_P,weight);

clearvars weight

%%
% Method 1
calcEnergyFlux(particle,angle);
% Method 2 (potentially wrong)
%calcEnergyFlux2(particle,angle,phi1,theta1,Bsurf,B_P);

% Test of mapping particles onto the surface
%mappingParticle(particle,angle,phi1,theta1,Bsurf,B_P);

