% Draw y=0 cuts for G8 and G28 flybys from GM 2D outputs.
%
% Hongyang Zhou, hyzhou@umich.edu

clear; clc
%% Parameters
% All the following sections need parameters here.
flyby = 'G8'; % default is G8

Rg = 2634000; %[m], radius of Ganymede

switch flyby
   case 'G8'
      disp('G8 flyby')
      filename = '~/Ganymede/MOP2018/runG8_PIC_1200s/GM/y=0_var_2_t00000520_n00229773.out'; % 2d GM outputs
      plotrange = [-3 3 -3 3];
   case 'G28'
      disp('G28 flyby')
      filename = '~/Ganymede/newPIC/run_G28_newPIC/x=0_var_1_n60000_263035.outs'; % 2d GM outputs
      plotrange = [-8 8 -8 8];
end

%% Postprocessing

[filehead,data] = read_data(filename,'verbose',false);
data = data.file1;


%% Plots

plot_data(data,filehead,'Pe Bx;Bz','plotrange',plotrange,...
  'plotmode','contf streamover','plotinterval',0.1);
axis equal; colorbar; caxis([0 2])

hold on
rectangle('Position',[-.5,-.5,1,1],'Curvature',[1,1]...
   ,'FaceColor',[.6 .6 .6]);
viscircles([0 0],1,'color',[.4 .2 1]);