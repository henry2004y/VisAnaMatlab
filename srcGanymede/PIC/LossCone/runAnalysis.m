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
% 11.Collect particles in sets of field lines.
%
% 11.Calculate flux from bulk velocity and pressure.
%
%
% Hongyang Zhou, hyzhou@umich.edu 10/10/2018

clear; clc; close all
%%

% Get pitch angles for all particles
[angle,B_P,particle,weight] = get_pitch_angle;

% Get particles inside the loss cone
[particle] = get_losscone(particle,angle,B_P);

% Get bulk velocity and pressure
p =  get_pressure(particle,weight);

%%
figure
scatter3(particle(1,:),particle(2,:),particle(3,:),'.');
axis equal