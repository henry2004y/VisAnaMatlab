% Comparison with Galileo PLS data along spacecraft trajectory.
% x=0 plane contour plot with field lines and quivers showing U.
% Note that the traj is not aligned with the plane, so it is actually a
% projection of velocity onto the plane.
%
% Hongyang Zhou, hyzhou@umich.edu 06/28/2018

clear; clc
%% Read Galileo trajectory data
flyby = 1;
flybyfile = strcat('Galileo_G',int2str(flyby),'_flyby_MAG.dat');
f = fullfile('~/Ganymede/GalileoData/galileomagdata',flybyfile);
[~,data] = read_log_data(f);

time = datetime(data(:,1:6));
xyz  = data(:,7:9);

%% Read simulation data
%filename='~/Documents/research/Ganymede/SteadyRun/run_G1_insulatingB1_5000/GM/x*';
filename='~/SWMF/SWMF/GMTSRUS/run_test/RESULTS/run_G1_test3/GM/x*outs';
npict = 21;
[filehead,data] = read_data(filename,'npict',npict);
data = data.file1;

% trajectory box output
filename='~/SWMF/SWMF/GM/BATSRUS/run_test/RESULTS/run_G1_test3/GM/box*outs';
[filehead_traj,data_traj] = read_data(filename,'npict',3);

data_traj = data_traj.file1;
x = data_traj.x(:,:,:,1);
y = data_traj.x(:,:,:,2);
z = data_traj.x(:,:,:,3);
ux_ = strcmpi('ux',filehead_traj.wnames);
uy_ = strcmpi('uy',filehead_traj.wnames);
uz_ = strcmpi('uz',filehead_traj.wnames);
ux = data_traj.w(:,:,:,ux_);
uy = data_traj.w(:,:,:,uy_);
uz = data_traj.w(:,:,:,uz_);

% From ndgrid to meshgrid format
x  = permute(x,[2 1 3]);
y  = permute(y,[2 1 3]);
z  = permute(z,[2 1 3]);
ux = permute(ux,[2 1 3]);
uy = permute(uy,[2 1 3]);
uz = permute(uz,[2 1 3]);

xyzPlot = xyz(1:150:end,:);
Usim = Inf(size(xyzPlot,1),3);
Usim(:,1) = interp3(x,y,z,ux,xyzPlot(:,1),xyzPlot(:,2),xyzPlot(:,3));
Usim(:,2) = interp3(x,y,z,uy,xyzPlot(:,1),xyzPlot(:,2),xyzPlot(:,3));
Usim(:,3) = interp3(x,y,z,uz,xyzPlot(:,1),xyzPlot(:,2),xyzPlot(:,3));

%% Plot
plotrange = [-4 4 -4 4];

plot_data(data,filehead,'uz by;bz','plotmode','contbar streamover',...
   'plotrange',plotrange,'plotinterval',0.05)
hold on
rectangle('Position',[-1.,-1.,2,2],'Curvature',[1,1]...
   ,'FaceColor',[.6 .6 .6]);

plot(xyz(:,2),xyz(:,3),'--r','LineWidth',2)
q = quiver(xyzPlot(:,2),xyzPlot(:,3),Usim(:,2),Usim(:,3));
q.Color = 'k';
q.LineWidth = 0.8;
