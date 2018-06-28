% Capturing magnetopause in 3D 
%
% Purpose: identify the boundary of status.
%
% This script is transformed into a function find_boundary_points.
%
% Hongyang Zhou, hyzhou@umich.edu 10/09/2017

clear; clc
%% Read in outputs
filename = '../run_G8_status_test/GM/3d_idl.outs'; % 3d GM outputs
[filehead,data] = read_data(filename,'verbose',false);


%% Find the positions of the status boundary
% The average of 0 and 2 is 1, so the boundary of 1 is not clear due to
% intepolation of the cut. If it is a truly 3D output, then there`s no such
% problem.
status = data.file1.w;
x = data.file1.x(:,:,:,1); 
y = data.file1.x(:,:,:,2);
z = data.file1.x(:,:,:,3);

%clearvars data filehead
%histogram(status)

% status=-3 are one the outer boundary or magnetopause or tail boundary;
% status=-2 is one point on the tail
% x1 = x(status==-2 );
% y1 = y(status==-2 );
% z1 = z(status==-2 );
% scatter3(x1,y1,z1);

status(status==-3) = 0;
status(status==-2) = 3;

% Pick the closed field line boundary points
x3_1 = x(status==3 & x<0);
y3_1 = y(status==3 & x<0);
z3_1 = z(status==3 & x<0);

x3_2 = x(status==3 & x>0);
y3_2 = y(status==3 & x>0);
z3_2 = z(status==3 & x>0);

k31   = boundary(x3_1,y3_1,z3_1,1);
k32   = boundary(x3_2,y3_2,z3_2,1);

bc31_ = unique(k31);
bc32_ = unique(k32);

%figure; hold on
%scatter3(x3_1(bc31_),y3_1(bc31_),z3_1(bc31_),'.');
%scatter3(x3_2(bc32_),y3_2(bc32_),z3_2(bc32_),'.');
%hold off;

% Pick only the dayside boundary points
x3_1 = x3_1(bc31_);
y3_1 = y3_1(bc31_);
z3_1 = z3_1(bc31_);

% Find the outer boundary points
mapindex_ = (x3_1.^2 + y3_1.^2 + z3_1.^2) > 1.81^2;
x3_1 = x3_1(mapindex_);
y3_1 = y3_1(mapindex_);
z3_1 = z3_1(mapindex_);

figure;
scatter3(x3_1,y3_1,z3_1,'.');

figure;
scatter(y3_1,z3_1,'.');

% Pick a boundary projection surface slightly larger than the last closed
% field lines.
figure;
scatter(1.05*y3_1,1.05*z3_1,'.');

return

figure; hold on
trisurf(k31,x3_1,y3_1,z3_1,'Facecolor','red','FaceAlpha',0.1,...
   'Edgecolor','k');
trisurf(k32,x3_2,y3_2,z3_2,'Facecolor','red','FaceAlpha',0.1,...
   'Edgecolor','k');
set(gca,'FontSize',16,'LineWidth',1.2)
hold off



%%
x1 = x(status==1);
y1 = y(status==1);
z1 = z(status==1);

x2 = x(status==2);
y2 = y(status==2);
z2 = z(status==2);

k1 = boundary(x1,y1,z1,1);
k2 = boundary(x2,y2,z2,1);
bc1_ = unique(k1);
bc2_ = unique(k2);

figure; hold on
scatter3(x1(bc1_),y1(bc1_),z1(bc1_),'.');
scatter3(x2(bc2_),y2(bc2_),z2(bc2_),'.');
scatter3(x3_1(bc31_),y3_1(bc31_),z3_1(bc31_),'.');
hold off; axis equal
xlim([-2 50]);ylim([-100 100]); zlim([-100 100]);
xlabel('x'); ylabel('y'); zlabel('z');
set(gca,'FontSize',16,'LineWidth',1.2)

% figure; hold on
% plot(x3_1(k31),z3_1(k31),'b'); 
% plot(x3_2(k32),z3_2(k32),'b');
% hold off; axis equal
