% E^\prime in Ohm`s law
%
% This should be collected into ReconnectionIdeaTest.m

clear; clc
%% Parameters
filenameGM = '~/Ganymede/newPIC/run_G8_newPIC/box_FTE_G8_1200s.outs';
ipict = 30;
s = 0.5; % compact boundary factor [0,1]
rThres = 1.6; % [Rg] radius threshold for finding dayside boundary pts
xThres = -1.5;

%% Find magnetopause, generate mesh

[x3bc,y3bc,z3bc] = ...
   find_boundary_points( filename,s,DoPlot,xThres,rThres,ipict );

[filehead,data] = read_data(filenameGM,'verbose',false,'npict',ipict);

% Fit the closed field line boundary with hypersurface

[fitresult,gof] = surface_fit(x3bc,y3bc,z3bc);

% Generate meshgrid from fitted surface
ymin = -1.; ymax = 1.; zmin = -0.5; zmax = 0.5;
dy = 1/30; dz = dy;
[yq,zq] = meshgrid(ymin:dy:ymax,zmin:dz:zmax);
xq = fitresult(yq,zq);

%% PC 
filenamePC='~/Ganymede/newPIC/run_G8_newPIC/3d_fluid_35.outs';
[filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);

data = data.file1;
x = data.x(:,:,:,1); 
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);
Bx = data.w(:,:,:,3);
By = data.w(:,:,:,4);
Bz = data.w(:,:,:,5);
Ex = data.w(:,:,:,6)*1e-3; % [mV/m]
Ey = data.w(:,:,:,7)*1e-3; % [mV/m]
Ez = data.w(:,:,:,8)*1e-3; % [mV/m]
Uex = data.w(:,:,:,9);
Uey = data.w(:,:,:,10);
Uez = data.w(:,:,:,11);
Uix = data.w(:,:,:,12);
Uiy = data.w(:,:,:,13);
Uiz = data.w(:,:,:,14);
Pe = data.w(:,:,:,15);
Pi = data.w(:,:,:,16);

% The original data is saved in ndgrid format. For streamline and
% isonormals functions, the input should be in meshgrid format.
x  = permute(x,[2 1 3]);
y  = permute(y,[2 1 3]);
z  = permute(z,[2 1 3]);
Bx = permute(Bx,[2 1 3]);
By = permute(By,[2 1 3]);
Bz = permute(Bz,[2 1 3]);
Ex = permute(Ex,[2 1 3]);
Ey = permute(Ey,[2 1 3]);
Ez = permute(Ez,[2 1 3]);
Uex = permute(Uex,[2 1 3]);
Uey = permute(Uey,[2 1 3]);
Uez = permute(Uez,[2 1 3]);
Uix = permute(Uix,[2 1 3]);
Uiy = permute(Uiy,[2 1 3]);
Uiz = permute(Uiz,[2 1 3]);
Pe  = permute(Pe,[2 1 3]);
Pi  = permute(Pi,[2 1 3]);

% Electron pressure gradient + electron intertia term [mV/m]
UeCrossB = -cross([Uex(:) Uey(:) Uez(:)],[Bx(:) By(:) Bz(:)])*1e-3;
Residual = [Ex(:) Ey(:) Ez(:)] - UeCrossB;

Residual = sqrt(sum(Residual.^2,2));

UeCrossBx = reshape(UeCrossB(:,1),size(Ex));
UeCrossBy = reshape(UeCrossB(:,2),size(Ex));
UeCrossBz = reshape(UeCrossB(:,3),size(Ex));
Residual = reshape(Residual,size(Ex));

[M,I] = max(Residual(:));
[Ix, Iy, Iz] = ind2sub(size(Residual),I);

% Maybe I need some smoothing: smooth3 for example
figure
a = slice(x,y,z,Residual,x(Ix,Iy,Iz),y(Ix,Iy,Iz),z(Ix,Iy,Iz)); 
colorbar; alpha(a,0.5)
a(1).LineStyle = 'none'; a(2).LineStyle = 'none'; a(3).LineStyle = 'none';
axis equal

hold on; scatter3(x(Ix,Iy,Iz),y(Ix,Iy,Iz),z(Ix,Iy,Iz),'r*')

%xRange = x(Ix,Iy,Iz)-0.1:0.1:x(Ix,Iy,Iz)+0.1;
xRange = x(Ix,Iy,Iz)-0.1;
yRange = y(Ix,Iy,Iz)-0.2:0.2:y(Ix,Iy,Iz)+1;
%zRange = z(Ix,Iy,Iz)-0.15:0.1:z(Ix,Iy,Iz)+0.15;
zRange = z(Ix,Iy,Iz);

[sx,sy,sz] = meshgrid(xRange, yRange, zRange);

hlines = streamline(x,y,z,Bx,By,Bz,sx,sy,sz);
%hlines = streamslice(x,y,z,Bx,By,Bz,sx(:),sy(:),sz(:));
set(hlines,'LineWidth',2,'Color','r');

hlines = streamline(x,y,z,-Bx,-By,-Bz,sx,sy,sz);
set(hlines,'LineWidth',2,'Color','r');

xlabel('x [R_G]'); ylabel('y [R_G]'); zlabel('z [R_G]');


%%
% The sliceomatic function is very good, but it can be improved:
% 1. colorbar range!
% sliceomatic(Pe)
% slicomatic(Residual)

%% Electron pressure gradient

e = 1.6022e-19; %[C]

RhoI = data.w(:,:,:,2)*1e6/14; %[/m^3]
PeXX = data.w(:,:,:,17); %[nPa]
PeYY = data.w(:,:,:,18);
PeZZ = data.w(:,:,:,19);
PeXY = data.w(:,:,:,20);
PeXZ = data.w(:,:,:,21);
PeYZ = data.w(:,:,:,22);

RhoI = permute(RhoI,[2 1 3]);
PeXX = permute(PeXX,[2 1 3]);
PeYY = permute(PeYY,[2 1 3]);
PeZZ = permute(PeZZ,[2 1 3]);
PeXY = permute(PeXY,[2 1 3]);
PeXZ = permute(PeXZ,[2 1 3]);
PeYZ = permute(PeYZ,[2 1 3]);

h = 1/32; % grid spacing
% Pressure gradient in x,y,z [mV/m]
GradPeX = 1e-12*(gradient(PeXX,h) + gradient(PeXY,h) + gradient(PeXZ,h)) ...
   ./ -(RhoI*e);
GradPeY = 1e-12*(gradient(PeXY,h) + gradient(PeYY,h) + gradient(PeYZ,h)) ...
   ./ -(RhoI*e);
GradPeZ = 1e-12*(gradient(PeXZ,h) + gradient(PeYZ,h) + gradient(PeZZ,h)) ...
   ./ -(RhoI*e);

%%
zplane = 77;

cutX = x(:,:,zplane);
cutY = y(:,:,zplane);

figure;
contourf(cutX,cutY,Ey(:,:,zplane));
colorbar; axis equal; hold on
scatter(cutX(:),cutY(:),'+')
title('Ey [mV/m]')

figure;
contourf(cutX,cutY,UeCrossBy(:,:,zplane));
colorbar; axis equal
title('U_e \times B [mV/m]')

figure;
contourf(cutX,cutY,GradPeY(:,:,zplane));
colorbar; axis equal; hold on
scatter(cutX(:),cutY(:),'+')
title('\nabla P_e/en_i [mV/m]')

figure;
contourf(cutX,cutY,Uex(:,:,zplane));
colorbar; axis equal; hold on
%scatter(cutX(:),cutY(:),'+')
quiver(cutX,cutY,Uex(:,:,zplane),Uey(:,:,zplane),2,'r');
title('U_{ex} [km/s]')

figure;
contourf(cutX,cutY,Uey(:,:,zplane));
colorbar; axis equal; hold on
%scatter(cutX(:),cutY(:),'+')
quiver(cutX,cutY,Uex(:,:,zplane),Uey(:,:,zplane),1,'r');
title('U_{ey} [km/s]')

figure;
contourf(cutX,cutY,Uiy(:,:,zplane));
colorbar; axis equal; hold on
%scatter(cutX(:),cutY(:),'+')
quiver(cutX,cutY,Uix(:,:,zplane),Uiy(:,:,zplane),1,'r');
title('U_{iy} [km/s]')
