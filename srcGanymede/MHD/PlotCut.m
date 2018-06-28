% Draw x=0 cuts for G8 and G28 flybys from GM 2D outputs.
%
% Hongyang Zhou, hyzhou@umich.edu

clear; clc
%% Parameters
% All the following sections need parameters here.
flyby = 'G28'; % default is G8

Rg = 2634000; %[m], radius of Ganymede

switch flyby
   case 'G8'
      disp('G8 flyby CPCP calculation')
      filename = '~/Ganymede/newPIC/run_G8_newPIC/x=0_var_1_n60000_247407.outs'; % 2d GM outputs
      plotrange = [-5 5 -8 8];
   case 'G28'
      disp('G28 flyby CPCP calculation')
      filename = '~/Ganymede/newPIC/run_G28_newPIC/x=0_var_1_n60000_263035.outs'; % 2d GM outputs
      plotrange = [-8 8 -8 8];
end

%% Postprocessing

[filehead,data] = read_data(filename,'verbose',false);
data = data.file1;

% non-integer status values need to be cleaned up for better visualization
status = data.w(:,:,14);
status(status>1 & status<2) = 2;
status(status<1) = 0;
data.w(:,:,14) = status;


%% Plots

plot_data(data,filehead,'status by;bz','plotrange',[-10 10 -10 10],...
  'plotmode','contf streamover','plotinterval',0.1);
axis equal
set(gca,'Xdir','reverse')

% func = 'ux by;bz';
% plotmode = 'contbar streamover';
% plot_data(data,filehead,func,'plotrange',plotrange,...
%    'plotmode',plotmode,'plotinterval',0.05);
% axis equal
% set(gca,'Xdir','reverse')