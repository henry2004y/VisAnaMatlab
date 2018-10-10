% Pitch angle calculation
%
% Hongyang Zhou, hyzhou@umich.edu 10/10/2018

clear; clc; %close all
%%

[xP,yP,zP,ux,uy,uz,weight] = get_particle;

[xF,yF,zF,Bx,By,Bz] = get_field;

[nP,angle,B_P,particle] = ...
   get_pitch_angle(xP,yP,zP,ux,uy,uz,xF,yF,zF,Bx,By,Bz);

[particle] = get_losscone(particle,angle,B_P,nP);