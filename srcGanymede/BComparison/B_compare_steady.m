% Script for SWMF Galileo Flyby Simulation Comparison
% Purpose: Plot quasi-steady state magnetic field comparison. Second norm
% is calculated to quantitatively estimate the difference.
%
% Prerequisite: Folder FileIO_SWMF must be included into Matlab path.
%
% Highlight feature:
% 1. finding the index for magnetic components automatically .
% 2. adding all kinds of interpolation methods available in Matlab. This
% part can work as a demo for using interpolation methods. Options:
%  (1) interp3
%  (2) griddata
%  (3) scatteredInterpolant
%  (4) griddedInterpolant
%
%
% Hongyang Zhou, hyzhou@umich.edu 08/08/2017
% Modified on 1/1/2018, 02/06/2018

%clear;clc
%% Parameters
flyby = 28;   % [1,2,7,8,28,29]
DoPlot = 1;  % Plot output
DoSave = 0;  % Save norm2 number
% Select the method for interpolation
InterpolationMethod = 1;

%% Read observation data
flybyfile = strcat('Galileo_G',int2str(flyby),'_flyby_MAG.dat');
f = fullfile('~/Ganymede/GalileoData/galileomagdata',flybyfile);
[~,data] = read_log_data(f);

time = datetime(data(:,1:6));
xyz  = data(:,7:9);
Bobs = data(:,10:12);
BobsStrength = sqrt(Bobs(:,1).^2+Bobs(:,2).^2+Bobs(:,3).^2);

%% Read/Plot simulation data  
%filename = strcat('~/Ganymede/newPIC/run_G8_newPIC/box_B_G8_1200s.outs');
%filename = strcat('~/Ganymede/newPIC/run_G28_newPIC/box_B_G28_1200s.outs');
%filename = strcat('~/Ganymede/newPIC/G2/box*.out');
filename = '~/Ganymede/MOP2018/runG*/GM/box*';
npict = 11; % Remember to change this for different runs!
[filehead,data] = read_data(filename,'npict',npict);

% Interpolate simulation data to observation data
data = data.file1;
nx = filehead.nx; nw = filehead.nw;

x = data.x(:,:,:,1);
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

w = data.w;
% Find the correct index of variables
bx_ = strcmpi('bx',filehead.wnames);
by_ = strcmpi('by',filehead.wnames);
bz_ = strcmpi('bz',filehead.wnames);
bx = data.w(:,:,:,bx_);
by = data.w(:,:,:,by_);
bz = data.w(:,:,:,bz_);

switch InterpolationMethod
   case {1,2,3}
      % From ndgrid to meshgrid format
      x  = permute(x,[2 1 3]);
      y  = permute(y,[2 1 3]);
      z  = permute(z,[2 1 3]);
      bx = permute(bx,[2 1 3]);
      by = permute(by,[2 1 3]);
      bz = permute(bz,[2 1 3]);
   case 4
      % use ndgrid
end


if InterpolationMethod==1
   % gridded interpolation, much faster than scattered interpolation
   Bsim = Inf(size(xyz,1),3);
   Bsim(:,1) = interp3(x,y,z,bx,xyz(:,1),xyz(:,2),xyz(:,3));
   Bsim(:,2) = interp3(x,y,z,by,xyz(:,1),xyz(:,2),xyz(:,3));
   Bsim(:,3) = interp3(x,y,z,bz,xyz(:,1),xyz(:,2),xyz(:,3));
elseif InterpolationMethod==2
   % Using griddata will only have interpolation
   % This is about 3 times faster than scatteredInterpolant
   x = data.x;
   x = reshape(x,[nx(1)*nx(2)*nx(3),3]);
   w = reshape(w,[nx(1)*nx(2)*nx(3),nw]);
   Bsim = Inf(size(xyz,1),3);
   
   Bsim(:,1) = griddata(x(:,1),x(:,2),x(:,3),w(:,bx_)...
      ,xyz(:,1),xyz(:,2),xyz(:,3),'linear');
   Bsim(:,2) = griddata(x(:,1),x(:,2),x(:,3),w(:,by_)...
      ,xyz(:,1),xyz(:,2),xyz(:,3),'linear');
   Bsim(:,3) = griddata(x(:,1),x(:,2),x(:,3),w(:,bz_)...
      ,xyz(:,1),xyz(:,2),xyz(:,3),'linear');
elseif InterpolationMethod==3
   % Using scatteredInterpolant will have both interpolation and 
   % extrapolation
   x = data.x;
   x = reshape(x,[nx(1)*nx(2)*nx(3),3]);
   w = reshape(w,[nx(1)*nx(2)*nx(3),nw]);
   Bsim = Inf(size(xyz,1),3);
   F1 = scatteredInterpolant(x(:,1),x(:,2),x(:,3),w(:,bx_));
   F1.Method = 'linear';
   F2 = scatteredInterpolant(x(:,1),x(:,2),x(:,3),w(:,by_));
   F2.Method = 'linear';
   F3 = scatteredInterpolant(x(:,1),x(:,2),x(:,3),w(:,bz_));
   F3.Method = 'linear';
   
   Bsim(:,1) = F1(xyz(:,1),xyz(:,2),xyz(:,3));
   Bsim(:,2) = F2(xyz(:,1),xyz(:,2),xyz(:,3));
   Bsim(:,3) = F3(xyz(:,1),xyz(:,2),xyz(:,3));
elseif InterpolationMethod==4
   % Using griddedInterpolant will have both interpolation and 
   % extrapolation

   Bsim = Inf(size(xyz,1),3);    
   F1 = griddedInterpolant(x,y,z,bx);
   F1.Method = 'linear';
   F2 = griddedInterpolant(x,y,z,by);
   F2.Method = 'linear';
   F3 = griddedInterpolant(x,y,z,bz);
   F3.Method = 'linear';
   
   Bsim(:,1) = F1(xyz(:,1),xyz(:,2),xyz(:,3));
   Bsim(:,2) = F2(xyz(:,1),xyz(:,2),xyz(:,3));
   Bsim(:,3) = F3(xyz(:,1),xyz(:,2),xyz(:,3));   
end

BsimStrength = sqrt(Bsim(:,1).^2+Bsim(:,2).^2+Bsim(:,3).^2);

%% Visualization
if DoPlot
   figure('Position', [100, 100, 1100, 800]);
   h = subplot(411); LW1 = 2.5; LW2 = 1.5; FS=16; Height = 0.18;
   plot(time, Bobs(:,1),'k',time,Bsim(:,1),'r','LineWidth',LW1); %Bx
   h.XTickLabel = []; 
   h.Position(4) = Height;
   ylabel('Bx [nT]');
   xlim([min(time) max(time)]);
   ylim([min(Bobs(:,1))-50 max(Bobs(:,1))+50 ]);
   %title({'(a)',strcat('Galileo G',int2str(flyby),' Flyby Magnetic field')})
   title(strcat('Galileo G',int2str(flyby),' Flyby Magnetic field'))
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(412);
   plot(time, Bobs(:,2),'k',time,Bsim(:,2),'r','LineWidth',LW1); %By
   h.XTickLabel = [];
   h.Position(4) = Height;   
   ylabel('By [nT]');
   xlim([min(time) max(time)]);
   ylim([min(Bobs(:,2))-50 max(Bobs(:,2))+50 ]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(413);
   plot(time, Bobs(:,3),'k',time, Bsim(:,3),'r','LineWidth',LW1); %Bz
   h.XTickLabel = [];
   h.Position(4) = Height;   
   ylabel('Bz [nT]');
   xlim([min(time) max(time)]);
   ylim([min(Bobs(:,3))-50 max(Bobs(:,3))+50 ]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
   h = subplot(414);   
   plot(time,BobsStrength,'k',time,BsimStrength,'r','LineWidth',LW1); %|B|
   h.Position(4) = Height;   
   ylabel('B [nT]');
   xlim([min(time) max(time)]);
   ylim([min(BobsStrength)-50 max(BobsStrength)+50 ]);
   set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
end

%% Calculate differences & save to output file
Bdiff = Bobs - Bsim;
% There`s a difference between 2 definitions of 2-norm.
%norm(Bdiff(~isnan(Bdiff)))
norm2 = mean(abs(Bdiff(~isnan(Bdiff))).^2)^(1/2);

if DoSave
   fileID = fopen('Bdiff.log','a');
   fprintf(fileID,'norm(Bdiff) %f\n',norm2);
   fclose(fileID);
end

% Find the index and value of the best parameter set
[val, idx] = min(norm2)