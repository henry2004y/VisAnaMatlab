% Galileo data processing
% Visualize the close encounter magnetic field and spacecraft trajectory. 
% It is useful to select PIC box location with this script.
% 
% Hongyang Zhou, hyzhou@umich.edu 01/16/2017
% Modified 06/28/2018
%-------------------------------------------------------------------------

clear; clc

%galileo_b_plot('Dir','../../GalileoData/galileomagdata','flyby','8','DoMovie',0,...
%   'DoPlotB',1,'DoSaveB',0);

galileo_b_plot('Dir','../../Galileo/','flyby','29');

%% G1
Start = [-0.090 -4.007 0.037];
End   = [2.055 3.531 1.203];
Center= 0.5*(Start+End);

Color = [1 0 0];
alpha = .6;
%plotcube(EdgeLength,Center-0.5*EdgeLength,alpha,Color)

plotcube(End-Start,Start,alpha,[0 1 0])

%% G2
Start = [-0.862 -4.696 0.581];
End   = [1.207 4.864 1.168];
Center= 0.5*(Start+End);


alpha = .6;
%plotcube(EdgeLength,Center-0.5*EdgeLength,alpha,Color)

plotcube(End-Start,Start,alpha,[0 1 0])

%% G7
Start = [0.915 4.919 1.749];
End   = [1.498 -5.764 1.697];
Center= 0.5*(Start+End);


Color = [1 0 0];
alpha = .6;

plotcube(End-Start,Start,alpha,Color)

%% G8
Start = [-1.575 -3.676 0.649];
End   = [-1.005 5.156 0.796];
Center= 0.5*(Start+End);


Color = [1 0 0];
alpha = .6;

plotcube(End-Start,Start,alpha,Color)

%% G28
Start = [-0.452 -7.700 -0.430];
End   = [-1.738 7.634 -0.330];
Center= 0.5*(Start+End);


Color = [0 1 0];
alpha = .6;


plotcube(End-Start,Start,alpha,Color)

%% G29
Start = [0.905 -7.084 1.628];
End   = [0.756 6.868 1.547];
Center= 0.5*(Start+End);


Color = [0 1 0];
alpha = .6;

plotcube(End-Start,Start,alpha,Color)