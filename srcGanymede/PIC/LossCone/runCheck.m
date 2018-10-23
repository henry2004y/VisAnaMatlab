% Script for checking the assumptions in particle flux calculation
%
% Do verification with uniform flux at the orginal locations!
%
% Hongyang Zhou, hyzhou@umich.edu 10/16/2018

clear; clc; %close all

%% Get pitch angles for all particles
[angle,Bx_P,By_P,Bz_P,B_P,particle,weight] = getParticleInfo;

%% Get particles inside the loss cone
[particle,angle,Bsurf,Bx_P,By_P,Bz_P,B_P,theta1,phi1] = ...
   getLossCone(particle,angle,Bx_P,By_P,Bz_P,B_P,weight);

clearvars weight Bx_P By_P Bz_P B_P theta1 phi1 angle

%%
calcEnergyFluxTest(particle);