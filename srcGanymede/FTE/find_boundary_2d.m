% Capturing magnetopause in 2D 
%
% Purpose: identify the boundary of status from 2D cuts.
%
% This script is transformed into a function find_boundary_points.
%
% Hongyang Zhou, hyzhou@umich.edu 10/08/2017

clear; clc
%% Read in outputs
filename = '../../newPIC/run_G8_newPIC/y*.outs'; % 2d GM outputs
[filehead,data] = read_data(filename,'verbose',false);


%% Plots for reference
func      = 'status';
plotmode  = 'contbar';
plotrange = [-5 15 -25 25];

plot_data(data.file1,filehead(1),func);

%% Find the positions of the status boundary
% The average of 0 and 2 is 1, so the boundary of 1 is not clear due to
% intepolation of the cut. If it is a truly 3D output, then there`s no such
% problem.
status = data.file1.w(:,:,21);
x = data.file1.x(:,:,1); 
z = data.file1.x(:,:,2);

% It seems that status==0.5 is the right on the magnetopause. below x=0
% status=2.5 is the boundary between closed field line and half open ones.
%histogram(status)

status(status==-1.5) = 0;
status(status==0.5)  = 1;
status(status==2.5)  = 2;

data.file1.w(:,:,21) = status;

%plot_data(data.file1,filehead,func,'plotrange',plotrange,'plotinterval',0.1);

% Pick the closed field line pts
x3_1 = x(status==3 & x<0);
z3_1 = z(status==3 & x<0);
x3_2 = x(status==3 & x>0);
z3_2 = z(status==3 & x>0);

x1 = x(status==1);
z1 = z(status==1);
x2 = x(status==2);
z2 = z(status==2);

k31   = boundary(x3_1,z3_1,1);
k32   = boundary(x3_2,z3_2,1);
k1   = boundary(x1,z1,1);
k2   = boundary(x2,z2,1);

figure; hold on
plot(x3_1(k31),z3_1(k31),'b'); 
plot(x3_2(k32),z3_2(k32),'b');
plot(x1(k1),z1(k1),'k'); 
plot(x2(k2),z2(k2),'r'); hold off; axis equal
