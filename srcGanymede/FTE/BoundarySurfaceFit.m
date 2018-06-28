% Magnetopause surface fit for upstream data analysis.
%
% Main features:
% # Find boundary points using status variable from BATSRUS;
% # Fit a high order polynomial surface from scattered pts;
% # Generate regular mesh from hypersurface model;
% # Include a rotation of the grid for G28 flyby.
%
% The content of this script is largely incorporated into other scripts.
% 
%
% Hongyang Zhou, hyzhou@umich.edu 02/16/2018

clear; clc
%% Find boundary points from steady state solution
flyby = 'G28';

switch flyby
   case 'G8'
      filename = '~/Ganymede/newPIC/run_G8_newPIC/3d_G8_steady.outs';
   case 'G28'
      filename = '~/Ganymede/newPIC/run_G28_newPIC/3d_G28_steady.out';
end

s = 0.5; % compact boundary factor in [0,1]

[x3bc,y3bc,z3bc] = find_boundary_points( filename,s );

%% Fit the closed field line boundary with hypersurface

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

% Plot fit with data.
figure(1); %hold on
h = plot( fitresult );
legend( h, 'poly5, x=x(y,z)', 'Location', 'NorthEast' );
% Label axes
xlabel('x [R_G]')
ylabel('y [R_G]')
zlabel('z [R_G]')
grid on; axis equal

xx = get(h, 'XData');
yy = get(h, 'YData');
zz = get(h, 'Zdata');
set(h, 'XData', zz, 'YData', xx, 'ZData', yy);

hold on;
scatter3(x3bc,y3bc,z3bc,20,'r','filled'); hold off
%axis tight
xlim([-2 0]); ylim([-1.5 1.5]); zlim([-0.6 0.8]);
set(gca,'FontSize',14,'LineWidth',1.2)

%% Generate mesh points from fitted surface

switch flyby
   case 'G8'
      % G8
      ymin = -1.1+1/15; ymax = 1.1-1/15; zmin = -0.54+1/15; zmax = 0.8-1/15;
   case 'G28'
      % G28
      ymin = -1.1+5/30; ymax = 1.1-5/30; zmin = -0.5; zmax = 0.75-3/30;
end

dy = 1/30; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

xq = fitresult(yq,zq);

if strcmp(flyby,'G28')
   % Try to rotate the ndgrid by 13 degrees
   % Define rotation angle
   theta = 13/180*pi;
   % Define rotation matrix
   rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
   
   temp = [yq(:),zq(:)]*rot.' ;
   Yrot = reshape(temp(:,1),[numel(yq),1]);
   Zrot = reshape(temp(:,2),[numel(yq),1]);
   
   Xrot = fitresult(Yrot,Zrot);
   
   xq = reshape(Xrot,size(yq));
   yq = reshape(Yrot,size(yq));
   zq = reshape(Zrot,size(yq));
end

%%
figure
surf(xq,yq,zq,'FaceAlpha',0.5,'FaceColor',...
   [0 0.450980392156863 0.741176470588235]);
axis equal
grid on; box on
xlabel('x [R_G]'); ylabel('y [R_G]'); zlabel('z [R_G]');
set(gca,'LineWidth',1.2,'FontSize',14);
view([-30.3 15.6]);


