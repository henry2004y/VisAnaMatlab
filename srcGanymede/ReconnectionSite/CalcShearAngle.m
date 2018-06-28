% Script for plotting magnetic shear angle at dayside magnetopause.
% Purpose: show that if there is a relation between shear angle and
% reconnection.
%
% From one talk at AGU 2017, I remembered one guy using multiple subplots
% at different locations showing the information in each square. This might
% be useful for me too.
%
% Based on this, I would say it is better to show the angle info with
% contour lines.
%
% Hongyang Zhou, hyzhou@umich.edu  11/02/2017

close all; clear; clc
%% Find boundary points from steady state solution
flyby = 'G8'; % default is G8

switch flyby
   case 'G8'
      % 3d GM outputs
      %filename = '~/Ganymede/newPIC/run_G8_newPIC/3d_G8_steady.outs';
      filename = '~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600s.outs'; %PC
      %filename='~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600to1200.outs';
      load('ReconnectionCount.mat')
      %load('ReconnectionCount_G8_600to1200.mat')
      xThres = -1.5;
      % dipole-direction unit vector
      unitDipole = [18 -51.82 716.8]/sqrt(18^2+51.82^2+716.8^2);
   case 'G28'
      %filename = '~/Ganymede/newPIC/run_G28_newPIC/3d_G28_steady.out';
      filename = '~/Ganymede/newPIC/G28_PIC_theta51/3d_fluid_600s.outs';
      load('ReconnectionCount_G28.mat')
      xThres = -1.75;
      % dipole-direction unit vector
      unitDipole = [19.26 -16.54 716.8]/sqrt(19.26^2+16.54^2+716.8^2);
end
      
% s = 0.5; % compact boundary factor [0,1]
% 
% [x3bc,y3bc,z3bc] = find_boundary_points( filename,s );

s = 1; % compact boundary factor [0,1]
[x3bc,y3bc,z3bc] = find_bz0_boundary( filename,s,xThres );


%% Fit the closed field line boundary with paraboloid

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

%% Generate mesh points from fitted surface
%ymin = -1.1+1/15; ymax = 1.1-1/15; zmin = -0.54+1/15; zmax = 0.8-1/15;
ymin = -1.2; ymax = 1.2; zmin = -0.8+1/16; zmax = 0.8;
dy = 1/32; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

xq = fitresult(yq,zq);

%% Transform into local coordicate system

% Calculate the normal direction to the fitted surface
[V, W] = differentiate(fitresult, yq, zq);

U = -ones(size(V));

% [U,V,W]
% figure(2);
% quiver3(xq,yq,zq,U,V,W,2,'color','r')

%% get the three local directions

% Initialize local vectors
dM = Inf(3,size(xq,1),size(xq,2));
dL = dM; dN = dM;

% Get the three local directions
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      dN(:,ix,iy) = [U(ix,iy) V(ix,iy) W(ix,iy)];
      dM(:,ix,iy) = cross(dN(:,ix,iy),unitDipole);
      dL(:,ix,iy) = cross(dM(:,ix,iy),dN(:,ix,iy));
      
      % Normalization
      dM(:,ix,iy) = dM(:,ix,iy) / norm(dM(:,ix,iy));
      dL(:,ix,iy) = dL(:,ix,iy) / norm(dL(:,ix,iy));
      dN(:,ix,iy) = dN(:,ix,iy) / norm(dN(:,ix,iy));
   end
end

%% Shear angle in the local coordinate system

% Offset by a distance of 1 cell (depends on grid resolution)
dn = dy;

switch flyby
   case 'G8'
      % box outputs
      %filename = '~/Ganymede/newPIC/run_G8_newPIC/box_FTE_G8_1200s.outs';
      filename = '~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600s.outs';
      %filename='~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600to1200.outs';
   case 'G28'
      % box outputs
      filename = '~/Ganymede/newPIC/G28_PIC_theta51/3d_fluid_600s.outs';
end
 
[filehead,data] = read_data(filename,'verbose',false,'npict',454);

data = data.file1;
x = data.x(:,:,:,1);
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

bx_ = strcmpi('bx',filehead.wnames);
by_ = strcmpi('by',filehead.wnames);
bz_ = strcmpi('bz',filehead.wnames);

bx = data.w(:,:,:,bx_);
by = data.w(:,:,:,by_);
bz = data.w(:,:,:,bz_);

% From ndgrid to meshgrid format
x  = permute(x,[2 1 3]);
y  = permute(y,[2 1 3]);
z  = permute(z,[2 1 3]);
bx = permute(bx,[2 1 3]);
by = permute(by,[2 1 3]);
bz = permute(bz,[2 1 3]);

% Get the starting surface points value
isurf = -5;
bxv= interp3(x, y, z, bx, xq + isurf*dn*squeeze(dN(1,:,:)),...
   yq + isurf*dn*squeeze(dN(2,:,:)),...
   zq + isurf*dn*squeeze(dN(3,:,:)));
byv= interp3(x, y, z, by, xq + isurf*dn*squeeze(dN(1,:,:)),...
   yq + isurf*dn*squeeze(dN(2,:,:)), ...
   zq + isurf*dn*squeeze(dN(3,:,:)));
bzv= interp3(x, y, z, bz, xq + isurf*dn*squeeze(dN(1,:,:)),...
   yq + isurf*dn*squeeze(dN(2,:,:)), ...
   zq + isurf*dn*squeeze(dN(3,:,:)));

bL = sum([bxv(:) byv(:) bzv(:)]'.*reshape(dL,[3,numel(xq)]));
bL = reshape(bL,size(xq));

offset = zeros(size(xq));
% Loop over the outer surfaces
for isurf = -4:10
   bxv= interp3(x, y, z, bx, xq + isurf*dn*squeeze(dN(1,:,:)),...
      yq + isurf*dn*squeeze(dN(2,:,:)),...
      zq + isurf*dn*squeeze(dN(3,:,:)));
   byv= interp3(x, y, z, by, xq + isurf*dn*squeeze(dN(1,:,:)),...
      yq + isurf*dn*squeeze(dN(2,:,:)), ...
      zq + isurf*dn*squeeze(dN(3,:,:)));
   bzv= interp3(x, y, z, bz, xq + isurf*dn*squeeze(dN(1,:,:)),...
      yq + isurf*dn*squeeze(dN(2,:,:)), ...
      zq + isurf*dn*squeeze(dN(3,:,:)));
   
   % Get bL on the outside surface
   bLOut = sum([bxv(:) byv(:) bzv(:)]'.*reshape(dL,[3,numel(xq)]));
   bLOut = reshape(bLOut,size(xq));
   
   % Assuming bL changes monotonically
   offset = offset + dn * (bL.*bLOut<0) .* (isurf-1+bL./(bL-bLOut));
   bL = bLOut;
end

% Get the true magnetopause position
xq = xq + offset.*squeeze(dN(1,:,:));
yq = yq + offset.*squeeze(dN(2,:,:));
zq = zq + offset.*squeeze(dN(3,:,:));

% Get B outside
bxv= interp3(x, y, z, bx, xq + 6*dn*squeeze(dN(1,:,:)),...
   yq + 3*dn*squeeze(dN(2,:,:)),...
   zq + 3*dn*squeeze(dN(3,:,:)));
byv= interp3(x, y, z, by, xq + 6*dn*squeeze(dN(1,:,:)),...
   yq + 3*dn*squeeze(dN(2,:,:)), ...
   zq + 3*dn*squeeze(dN(3,:,:)));
bzv= interp3(x, y, z, bz, xq + 6*dn*squeeze(dN(1,:,:)),...
   yq + 3*dn*squeeze(dN(2,:,:)), ...
   zq + 3*dn*squeeze(dN(3,:,:)));

bLOut = Inf(size(xq)); bMOut = bLOut; %bN = bL;

% This could potentially be improved!
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      bLOut(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
      bMOut(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
      %bN(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dN(:,ix,iy);
   end
end


% Get B inside
bxv= interp3(x, y, z, bx, xq - 4*dn*squeeze(dN(1,:,:)),...
   yq - 3*dn*squeeze(dN(2,:,:)),...
   zq - 3*dn*squeeze(dN(3,:,:)));
byv= interp3(x, y, z, by, xq - 4*dn*squeeze(dN(1,:,:)),...
   yq - 3*dn*squeeze(dN(2,:,:)), ...
   zq - 3*dn*squeeze(dN(3,:,:)));
bzv= interp3(x, y, z, bz, xq - 4*dn*squeeze(dN(1,:,:)),...
   yq - 3*dn*squeeze(dN(2,:,:)), ...
   zq - 3*dn*squeeze(dN(3,:,:)));

bLIn = Inf(size(xq)); bMIn = bLIn; %bN = bL;

% This could potentially be improved!
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      bLIn(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
      bMIn(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
      %bN(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dN(:,ix,iy);
   end
end

% Calculate the shear angle
% equivalent to: acosd( dot(u,v)/(norm(u)*norm(v)) );
theta = acosd( (bLIn.*bLOut + bMIn.*bMOut) ./ ...
   sqrt((bLIn.^2+bMIn.^2).*(bLOut.^2+bMOut.^2)) );

figure('position', [10, 10, 800, 400])
contourf(yq,zq,theta,'Linestyle','none'); 
c = colorbar; colormap jet; caxis([0 180])
c.Label.String = 'degree';
axis equal
% hold on
% quiver(yq,zq,-sind(theta),cosd(theta)); hold off
% xlabel('y [R_G]'); ylabel('z [R_G]');
set(gca,'Xdir','reverse','FontSize',16,'LineWidth',1.2);

figure('position', [10, 10, 800, 400]);
[C,h] = contour(yq,zq,theta);
clabel(C,h,'Color','k');
h.ContourZLevel = 200;
switch flyby
   case 'G8'
      h.LevelList = [175 170 165 160 155 125 135 125]; % G8
   case 'G28'
      h.LevelList = [110 100 90 80 70 60 50 40]; % G8
end

colormap(flipud(gray));
%colormap gray


Probability = ReconnectionCount ./ 600;
hold on; surf(yq,zq,Probability,'LineStyle','none'); 
switch flyby
   case 'G8'
      caxis([0 0.35]) % G8
   case 'G28'
      caxis([0 0.16])  % G28
end

c = colorbar; c.Label.String = 'Normalized probability';
%c.Limits = [0 0.35]; % G8

index = ReconnectionCount>=20;
%scatter(yq(index),zq(index),'.','MarkerFaceColor',[0.5 0.5 0.5],...
%   'MarkerEdgeColor',[0.5 0.5 0.5])
scatter(yq(index),zq(index),'.')
xlabel('y [R_G]'); ylabel('z [R_G]')
axis equal
set(gca,'Xdir','reverse','FontSize',16,'LineWidth',1.2);


return
%% Shear angle in the local coordinate system as a movie

Coefin  = 0.9; % shrink factor
Coefout = 1.1; % expansion factor

switch flyby
   case 'G8'
      % box outputs
      filename = '~/Ganymede/newPIC/run_G8_newPIC/box_FTE_G8_1200s.outs'; 
   case 'G28'
      % box outputs
      filename = '~/Ganymede/newPIC/run_G28_newPIC/box_FTE_G28_1200s.outs';
end
 
v = VideoWriter('angle_test.avi');
v.FrameRate = 10;
v.open

fig = figure('position', [10, 10, 800, 600]); 
set(fig,'nextplot','replacechildren');

for ipict = 1:30
   [filehead,data] = read_data(filename,'verbose',false,'npict',ipict);
   
   x = data.file1.x(:,:,:,1);
   y = data.file1.x(:,:,:,2);
   z = data.file1.x(:,:,:,3);
   
   bx = data.file1.w(:,:,:,5);
   by = data.file1.w(:,:,:,6);
   bz = data.file1.w(:,:,:,7);
   
   % From ndgrid to meshgrid format
   x  = permute(x,[2 1 3]);
   y  = permute(y,[2 1 3]);
   z  = permute(z,[2 1 3]);
   bx = permute(bx,[2 1 3]);
   by = permute(by,[2 1 3]);
   bz = permute(bz,[2 1 3]);
   
   % Calculate the inner direction tangential to the boundary surface
   Coef = Coefin;
   bxv= interp3(x, y, z, bx, Coef*xq, Coef*yq, Coef*zq);
   byv= interp3(x, y, z, by, Coef*xq, Coef*yq, Coef*zq);
   bzv= interp3(x, y, z, bz, Coef*xq, Coef*yq, Coef*zq);
   
   % Transform vectors into local coordinate system
   b1in = Inf(size(xq)); b2in = b1in;
   
   % This could potentially be improved!
   for ix=1:size(xq,1)
      for iy=1:size(xq,2)
         b1in(ix,iy)  = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*d1(:,ix,iy);
         b2in(ix,iy)  = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*d2(:,ix,iy);
      end
   end
   
   % Calculate the outer direction tangential to the boundary surface
   Coef = Coefout;
   bxv= interp3(x, y, z, bx, Coef*xq, Coef*yq, Coef*zq);
   byv= interp3(x, y, z, by, Coef*xq, Coef*yq, Coef*zq);
   bzv= interp3(x, y, z, bz, Coef*xq, Coef*yq, Coef*zq);
   
   
   % Transform vectors into local coordinate system
   b1out = Inf(size(xq)); b2out = b1out;
   
   % This could potentially be improved!
   for ix=1:size(xq,1)
      for iy=1:size(xq,2)
         b1out(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*d1(:,ix,iy);
         b2out(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*d2(:,ix,iy);
      end
   end
   
   % Calculate the shear angle
   % equivalent to: acosd( dot(u,v)/(norm(u)*norm(v)) );
   theta = acosd( (b1in.*b1out + b2in.*b2out) ./ ...
      sqrt((b1in.^2+b2in.^2).*(b1out.^2+b2out.^2)) );
   

   contourf(yq,zq,theta,'Linestyle','none');
   c = colorbar; colormap jet; 
   c.Label.String = 'degree';
   axis equal; caxis manual; caxis([0 180])
   hold on
   quiver(yq,zq,-sind(theta),cosd(theta)); hold off
   xlabel('y [R_G]'); ylabel('z [R_G]');
   set(gca,'Xdir','reverse','FontSize',16,'LineWidth',1.2);
   
   
   frame = getframe(gcf);
   clf;
   writeVideo(v,frame);
end

v.close