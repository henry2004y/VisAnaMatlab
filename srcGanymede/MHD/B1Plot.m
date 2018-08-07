% Script for showing the differences in B1 between different resistivity
% setup, MOP2018
%
% Hongyang Zhou, hyzhou@umich.edu, 07/05/2018

clear; clc; close all
%%
%filename1 = '~/Documents/research/Ganymede/SteadyRun/Run34_floatB1_OldGrid_CoarseAxis/GM/y*';
%filename2 = '~/Document/research/Ganymede/SteadyRun/MOP2018/runG8_steady/GM/y*';

filename = ['~/Documents/research/Ganymede/SteadyRun/'...
   'Run34_floatB1_OldGrid_CoarseAxis/GM/y* '...
   '~/Documents/research/Ganymede/SteadyRun/MOP2018/runG8_steady/GM/y*'];

plotrange = [-2 2 -2 2];
func = 'b1z b1x;b1z';
plotmode = 'contbar streamover';
plotinterval = 0.1;

[filehead,data] = read_data(filename,'verbose',false,'npict',51);


%%
plot_data(data.file1,filehead(1),func,...
   'plotrange',plotrange,...
   'plotmode',plotmode,...
   'plotinterval',plotinterval)

caxis([-500 40]); %colormap jet

hold on
rectangle('Position',[-.5,-.5,1,1],'Curvature',[1,1]...
   ,'FaceColor',[.6 .6 .6]);
%viscircles([0 0],1,'color',[.4 .2 1]);
hold off


plot_data(data.file2,filehead(2),func,...
   'plotrange',plotrange,...
   'plotmode',plotmode,...
   'plotinterval',plotinterval)

caxis([-500 40]); %colormap jet

hold on
rectangle('Position',[-.5,-.5,1,1],'Curvature',[1,1]...
   ,'FaceColor',[.6 .6 .6]);
%viscircles([0 0],1,'color',[.4 .2 1]);
hold off