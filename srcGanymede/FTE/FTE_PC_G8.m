% Magnetopause contour plots animation from 3D PC outputs for G8.
%
% Repeat the same procedure of analyzing GM outputs as in G8FTE_LMN.m
% to PC outputs.
% From PC outputs, we can get electron velocity information.
% My impression from doing this is that Ve is not clear in showing the
% reconnection at the upstream magnetopause.
% This script needs to be improved for further usage.
%
% Hongyang Zhou, hyzhou@umich.edu
%
% modified 06/29/2018, version 1.2

clear; clc
%% Parameters
% 3D GM outputs
filename = '~/Ganymede/newPIC/run_G8_newPIC/3d_G8_steady.outs';
s = 0.8; % compact boundary factor [0,1]

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
VA  = 230; %[km/s], for normalization

%% Find boundary points from steady state solution

[x3bc,y3bc,z3bc] = find_boundary_points( filename,s );

%% Fit the closed field line boundary with hypersurface

[fitresult,gof] = surface_fit(x3bc,y3bc,z3bc);

%% Generate mesh points from fitted surface and calculate LMN directions
ymin = -1.2; ymax = 1.2; zmin = -0.8; zmax = 0.8;
dy = 1/32; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

xq = fitresult(yq,zq);

% Calculate the normal direction to the fitted surface
[V, W] = differentiate(fitresult, yq, zq);
U = -ones(size(V));

% get the three local directions
% dipole-direction unit vector
unitDipole = [18 -51.82 716.8]/sqrt(18^2+51.82^2+716.8^2);
% Initialize local vectors: d1-> M d2->L d3-> N
dL = Inf(3,size(xq,1),size(xq,2));
dM = dL; dN = dL;

% This part could potentially be optimized!
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      dN(:,ix,iy) = [U(ix,iy) V(ix,iy) W(ix,iy)];
      dM(:,ix,iy) = cross(dN(:,ix,iy),unitDipole);
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

% Create video
v = VideoWriter(Vname);
v.FrameRate = vFrameRate;
v.open

% create new figure with specified size
hfig = figure(4);
set(hfig,'position', [10, 10, 800, 600]) 
colormap(jet);

for ipict=1:10%npict
   fprintf('ipict=%d\n',ipict)
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
 
   
   % Find the relative position where BL goes to 0   
   bxvOut= interp3(x, y, z, bx, Coef*xq, Coef*yq, Coef*zq);
   byvOut= interp3(x, y, z, by, Coef*xq, Coef*yq, Coef*zq); 
   bzvOut= interp3(x, y, z, bz, Coef*xq, Coef*yq, Coef*zq);
   bxv= interp3(x, y, z, bx, xq, yq, zq);
   byv= interp3(x, y, z, by, xq, yq, zq); 
   bzv= interp3(x, y, z, bz, xq, yq, zq);   
   bL = Inf(size(xq));

   for ix=1:size(xq,1)
      for iy=1:size(xq,2)        
         bL(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
      end
   end  
   
   %
   nev= interp3(x, y, z, ne, Coef*xq, Coef*yq, Coef*zq);
   niv= interp3(x, y, z, ni, Coef*xq, Coef*yq, Coef*zq);
   bxv= interp3(x, y, z, bx, Coef*xq, Coef*yq, Coef*zq);
   byv= interp3(x, y, z, by, Coef*xq, Coef*yq, Coef*zq); 
   bzv= interp3(x, y, z, bz, Coef*xq, Coef*yq, Coef*zq);
   uixv= interp3(x, y, z, uix, Coef*xq, Coef*yq, Coef*zq);
   uiyv= interp3(x, y, z, uiy, Coef*xq, Coef*yq, Coef*zq); 
   uizv= interp3(x, y, z, uiz, Coef*xq, Coef*yq, Coef*zq);
   uexv= interp3(x, y, z, uex, Coef*xq, Coef*yq, Coef*zq);
   ueyv= interp3(x, y, z, uey, Coef*xq, Coef*yq, Coef*zq); 
   uezv= interp3(x, y, z, uez, Coef*xq, Coef*yq, Coef*zq);   
   piv= interp3(x, y, z, pi, Coef*xq, Coef*yq, Coef*zq);
   pev= interp3(x, y, z, pe, Coef*xq, Coef*yq, Coef*zq);
   

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
   % Note: This normalization is wrong! I shouldn't use different Alfven
   % velocity for different cells!
   
   % Current density
   J = e*(niv(:).*[uixv(:) uiyv(:) uizv(:)]./mi - ...
      nev(:).*[uexv(:) ueyv(:) uezv(:)])*1e9; %[\mu A/m^2]
   J = sqrt(sum(J.^2,2));
   J = reshape(J,size(xq));
   
   subplot_tight(4,3,1);
   contourf(yq,zq,uiL./VA,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-1 1])
   %caxis([-200 200])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('U_{iL} [km/s]')
   
   subplot_tight(4,3,2);
   contourf(yq,zq,uiM./VA,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-1.5 0.5])
   %caxis([-200 100])
   %xlabel('y [R_G]');
   ylabel('z [R_G]');
   title('U_{iM} [km/s]')
       
   subplot_tight(4,3,3);
   contourf(yq,zq,uiN./VA,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-0.5 0.5])
   %caxis([-60 100])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('U_{iN} [km/s]')
   
   % ueL
   subplot_tight(4,3,4);
   contourf(yq,zq,ueL./VA,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-4 4])
   %caxis([-600 600])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('U_{eL} [km/s]')
   
   % ueM
   subplot_tight(4,3,5);
   contourf(yq,zq,ueM./VA,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0 6])
   %caxis([0 1300])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('U_{eM} [km/s]')
   
   % ueN
   subplot_tight(4,3,6);
   contourf(yq,zq,ueN./VA,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-1 1])
   %caxis([-400 300])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('U_{eN} [km/s]')
   
   % BL (positive pointing upstream)
   subplot_tight(4,3,7);
   contourf(yq,zq,bL,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   %caxis([0 180])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('B_L [nT]')
   
   % BM (positive pointing in roughly y direction)
   subplot_tight(4,3,8);
   contourf(yq,zq,bM,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-80 40])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('B_M [nT]')
   
   % BN
   subplot_tight(4,3,9);
   contourf(yq,zq,bN,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-100 100])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('B_N [nT]')
   
   subplot_tight(4,3,10);
   contourf(yq,zq,piv,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([2 11])
   %xlabel('y [R_G]'); ylabel('z [R_G]');
   title('P_{i} [nPa]')
   
   subplot_tight(4,3,11);
   contourf(yq,zq,pev,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0.1 3.5])
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('P_{e} [nPa]')
   
   subplot_tight(4,3,12);
   contourf(yq,zq,J,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0 0.2])
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
