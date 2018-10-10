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
%
% 11.flux? pressure, bulk velocity
%
%
% Hongyang Zhou, hyzhou@umich.edu 10/10/2018

clear; clc; close all
%%

[xP,yP,zP,ux,uy,uz,weight] = get_particle;

[xF,yF,zF,Bx,By,Bz] = get_field;

[nP,angle,B_P,particle] = ...
   get_pitch_angle(xP,yP,zP,ux,uy,uz,xF,yF,zF,Bx,By,Bz);

[particle] = get_losscone(particle,angle,B_P,nP);

%%
figure
scatter3(particle(1,:),particle(2,:),particle(3,:),'.');
axis equal