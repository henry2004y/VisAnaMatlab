% Magnetopause contour plots movie from 3D PC for G28.
%
% The magnetopause is found by identifying the nearest location where BL
% changes sign.
% Approximate boundary points is found by tracking Bz=0.
%
% The orientation of LMN needs to be checked!!!
%
% Hongyang Zhou, hyzhou@umich.edu
%
% modified 06/29/2018, version 1.2

clear; clc
%% Parameters
% 3D PC outputs
filename = '~/Ganymede/newPIC/G28_PIC_theta51/3d_fluid_600s.outs';
s = 1.0; % compact boundary factor [0,1]
xThres = -1.75;

% GM upstream box outputs
filenamePC = '~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600s.outs';
% Output movie
Vname = 'test.avi';
vFrameRate = 10;

% Criteria for surface contour and FTE identification
Coef         = 1.08; % expansion factor from original surface fit 
threshold_pe = 2.1;
threshold_j  = 0.52; 

%% Physical Parameters
mu0 = 4*pi*1e-7;    %[H/m]
me  = 9.1094e-31;   %[kg] 
mp  = 1.6726*1e-27; %[kg]
mi  = 14; % average ion mass [amu]
e   = 1.6022e-19; %[C]
VA  = 450; %[km/s] Estimation of Alfven velocity for nomalization

%% Find boundary points from steady state solution

%[x3bc,y3bc,z3bc] = find_boundary_points( filename,s ); 
[x3bc,y3bc,z3bc] = find_bz0_boundary( filename,s,xThres );

%% Fit the closed field line boundary

% Set up fittype and options.
ft = fittype( 'poly55' );
%ft = fittype( 'thinplateinterp' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

%% Generate mesh points from fitted surface and calculate LMN directions
ymin = -1.15; ymax = 1.15; zmin = -0.75; zmax = 0.75;
dy = 1/32; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

% Try to rotate the ndgrid by 15 degrees
% Define rotation angle
theta = 15/180*pi; 
% Define rotation matrix
rot = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 

temp = [yq(:),zq(:)]*rot.' ;
Yrot = reshape(temp(:,1),[numel(yq),1]);
Zrot = reshape(temp(:,2),[numel(yq),1]);

Xrot = fitresult(Yrot,Zrot);

xq = reshape(Xrot,size(yq));
yq = reshape(Yrot,size(yq));
zq = reshape(Zrot,size(yq));

% xq = fitresult(yq,zq);
xq(xq>-1.13) = nan;

% Calculate the normal direction to the fitted surface
[V, W] = differentiate(fitresult, yq, zq);
U = -ones(size(V));

% get the three local directions
% dipole-direction unit vector
unitDipole = [19.26 -16.54 716.8]/sqrt(19.26^2+16.54^2+716.8^2);
% Upstream B unit vector
%unitUpstreamB = [-7 78 -76]/sqrt(7^2+78^2+76^2);

unitL = [0 -sind(15) cosd(15)];

% Initialize local vectors: d1-> M d2->L d3-> N
dL = Inf(3,size(xq,1),size(xq,2));
dM = dL; dN = dL;

% This part could potentially be optimized!
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      dN(:,ix,iy) = [U(ix,iy) V(ix,iy) W(ix,iy)];
      %dM(:,ix,iy) = cross(dN(:,ix,iy),unitDipole);
      dM(:,ix,iy) = cross(dN(:,ix,iy),unitL);
      dL(:,ix,iy) = cross(dM(:,ix,iy),dN(:,ix,iy));
      
      % Normalization
      dL(:,ix,iy) = dL(:,ix,iy) / norm(dL(:,ix,iy));
      dM(:,ix,iy) = dM(:,ix,iy) / norm(dM(:,ix,iy));
      dN(:,ix,iy) = dN(:,ix,iy) / norm(dN(:,ix,iy));
   end
end

%% Movie from PIC outputs

[~,~,fileinfo] = read_data(filenamePC,'verbose',false);
npict = fileinfo.npictinfiles;

% Offset by a distance of 1 cell (dependent on grid resolution)
dn = dy;

% Create video
v = VideoWriter(Vname);
v.FrameRate = vFrameRate;
v.open

% create new figure with specified size
hfig = figure(4);
set(hfig,'position', [10, 10, 800, 520]) 
colormap(jet);

for ipict = 1:npict
   [filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);
   
   data = data.file1;
   x = data.x(:,:,:,1);
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);

   ne = data.w(:,:,:,1)*1e6;    % [/m^3]
   ni = data.w(:,:,:,2)*1e6;    % [/m^3]
   bx   = data.w(:,:,:,3);      % [nT]
   by   = data.w(:,:,:,4);
   bz   = data.w(:,:,:,5);
   uex  = data.w(:,:,:,9);      % [km/s]
   uey  = data.w(:,:,:,10);
   uez  = data.w(:,:,:,11);
   uix  = data.w(:,:,:,12);
   uiy  = data.w(:,:,:,13);
   uiz  = data.w(:,:,:,14);
   pe   = data.w(:,:,:,15);     % [nPa]
   pi   = data.w(:,:,:,16);
   
   % From ndgrid to meshgrid format
   ne = permute(ne,[2 1 3]);
   ni = permute(ni,[2 1 3]);
   x  = permute(x,[2 1 3]);
   y  = permute(y,[2 1 3]);
   z  = permute(z,[2 1 3]);
   bx = permute(bx,[2 1 3]);
   by = permute(by,[2 1 3]);
   bz = permute(bz,[2 1 3]);   
   uex = permute(uex,[2 1 3]);
   uey = permute(uey,[2 1 3]);
   uez = permute(uez,[2 1 3]);
   uix = permute(uix,[2 1 3]);
   uiy = permute(uiy,[2 1 3]);
   uiz = permute(uiz,[2 1 3]);   
   pe = permute(pe,[2 1 3]);
   pi = permute(pi,[2 1 3]);
 
   
   % Get the starting surface points value
%    bxv= interp3(x, y, z, bx, xq, yq, zq);
%    byv= interp3(x, y, z, by, xq, yq, zq);
%    bzv= interp3(x, y, z, bz, xq, yq, zq);
   
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
   % Loop over the inner/outer surfaces
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
      %mesh(offset)
   end
   
   offset(offset==0) = nan; % Remove the unecessary points
   %offset = 0;
   xqNew = xq + offset.*squeeze(dN(1,:,:));
   yqNew = yq + offset.*squeeze(dN(2,:,:));
   zqNew = zq + offset.*squeeze(dN(3,:,:)); 
   
   nev= interp3(x, y, z, ne, xqNew, yqNew, zqNew);
   niv= interp3(x, y, z, ni, xqNew, yqNew, zqNew);
   bxv= interp3(x, y, z, bx, xqNew, yqNew, zqNew);
   byv= interp3(x, y, z, by, xqNew, yqNew, zqNew); 
   bzv= interp3(x, y, z, bz, xqNew, yqNew, zqNew);
   uixv= interp3(x, y, z, uix, xqNew, yqNew, zqNew);
   uiyv= interp3(x, y, z, uiy, xqNew, yqNew, zqNew); 
   uizv= interp3(x, y, z, uiz, xqNew, yqNew, zqNew);
   uexv= interp3(x, y, z, uex, xqNew, yqNew, zqNew);
   ueyv= interp3(x, y, z, uey, xqNew, yqNew, zqNew); 
   uezv= interp3(x, y, z, uez, xqNew, yqNew, zqNew);   
   piv= interp3(x, y, z, pi, xqNew, yqNew, zqNew);
   pev= interp3(x, y, z, pe, xqNew, yqNew, zqNew);
   
   % Transform vectors into local coordinate system
   uiL = Inf(size(xq)); uiM = uiL; uiN = uiL;
   ueL = Inf(size(xq)); ueM = ueL; ueN = ueL;
   bL = Inf(size(xq)); bM = bL; bN = bL;
   
   % This could potentially be improved!
   for ix=1:size(xq,1)
      for iy=1:size(xq,2)
         uiL(ix,iy) = [uixv(ix,iy) uiyv(ix,iy) uizv(ix,iy)]*dL(:,ix,iy);
         uiM(ix,iy) = [uixv(ix,iy) uiyv(ix,iy) uizv(ix,iy)]*dM(:,ix,iy);
         uiN(ix,iy) = [uixv(ix,iy) uiyv(ix,iy) uizv(ix,iy)]*dN(:,ix,iy);
         ueL(ix,iy) = [uexv(ix,iy) ueyv(ix,iy) uezv(ix,iy)]*dL(:,ix,iy);
         ueM(ix,iy) = [uexv(ix,iy) ueyv(ix,iy) uezv(ix,iy)]*dM(:,ix,iy);
         ueN(ix,iy) = [uexv(ix,iy) ueyv(ix,iy) uezv(ix,iy)]*dN(:,ix,iy);
         bL(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
         bM(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
         bN(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dN(:,ix,iy);
      end
   end
   
   % Alfven velocity
   %VA = sqrt((bL.^2 + bM.^2 + bN.^2) ./ (mu0*(niv*mp+nev*me)))*1e-12;  

   % Current density
   J = e*(niv(:).*[uixv(:) uiyv(:) uizv(:)]./mi - ...
      nev(:).*[uexv(:) ueyv(:) uezv(:)])*1e9; %[\mu A/m^2]
   J = sqrt(sum(J.^2,2));
   J = reshape(J,size(xq));
   
   subplot_tight(4,3,1);
   contourf(yqNew,zqNew,uiL./VA,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-0.5 0.5])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('U_{iL} [km/s]')
   
   subplot_tight(4,3,2);
   contourf(yqNew,zqNew,uiM./VA,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-0.6 0.1])
   %caxis([-200 100])
   %xlabel('y [R_G]');
   ylabel('z [R_G]');
   title('U_{iM} [km/s]')
       
   subplot_tight(4,3,3);
   contourf(yqNew,zqNew,uiN./VA,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-0.3 0.1])
   %caxis([-60 100])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('U_{iN} [km/s]')
   
   % ueL
   subplot_tight(4,3,4);
   contourf(yqNew,zqNew,ueL./VA,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-1.5 1.5])
   %caxis([-600 600])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('U_{eL} [km/s]')
   
   % ueM
   subplot_tight(4,3,5);
   contourf(yqNew,zqNew,ueM./VA,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-1 4])
   %caxis([0 1300])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('U_{eM} [km/s]')
   
   % ueN
   subplot_tight(4,3,6);
   contourf(yqNew,zqNew,ueN./VA,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-0.6 1])
   %caxis([-400 300])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('U_{eN} [km/s]')
   
   % BL (positive pointing upstream)
   subplot_tight(4,3,7);
   contourf(yqNew,zqNew,bL,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-10 10])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('B_L [nT]')
   
   % BM (positive pointing in roughly y direction)
   subplot_tight(4,3,8);
   contourf(yqNew,zqNew,bM,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([40 180])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('B_M [nT]')
   
   % BN
   subplot_tight(4,3,9);
   contourf(yqNew,zqNew,bN,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-60 60])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('B_N [nT]')
   
   % Pi
   subplot_tight(4,3,10);
   contourf(yqNew,zqNew,piv,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0.5 6])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('P_{i} [nPa]')
   
   % Pe
   subplot_tight(4,3,11);
   contourf(yqNew,zqNew,pev,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0.2 2.5])
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('P_{e} [nPa]')
   
   % J
   subplot_tight(4,3,12);
   contourf(yqNew,zqNew,J,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0 0.20])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('J [\mu A/m^2]')
    
   
   %
   dim = [0.2 0.01 0.05 0.02];
   str = sprintf('it=%d, time=%.1fs',filehead.it,filehead.time);
   a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
      'FontWeight','bold','EdgeColor','none');
   
   frame = getframe(gcf);
   
   set(gca,'nextplot','replacechildren');
   
   clf;
   writeVideo(v,frame);   
   
end

v.close
close(4)
