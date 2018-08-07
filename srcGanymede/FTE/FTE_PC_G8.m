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
IsGatheredFile = false; % Single or multiple input files
% 3D GM outputs
%filename = '~/Ganymede/newPIC/run_G8_newPIC/3d_G8_steady.outs';
%filename = '~/Ganymede/MOP2018/runG8_PIC_1200s/3d_t=0.out';
filename = '~/Ganymede/MOP2018/runG8_PIC_1200s/3d_t=280.out';
s = 0.8; % compact boundary factor [0,1]

% PC upstream box outputs
PCdir = '~/Ganymede/MOP2018/runG8_PIC_1200s/PC';
%filenamePC = '~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600s.outs';
%filenamePC = '3d_fluid_600s.outs';
% Output movie
Vname = '~/Ganymede/MOP2018/2DMagnetopause_G8.avi';
vFrameRate = 10;

% Criteria for surface contour and FTE identification
Coef         = 1.00; % expansion factor from original surface fit 
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
ymin = -1.2; ymax = 1.2; zmin = -0.6; zmax = 0.6;
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
if IsGatheredFile
   % single input file case
   [filehead,~,fileinfo] = read_data(filename,'verbose',false);
   npict = fileinfo.npictinfiles; % # of snapshot in the file
else
   % multiple input file case
   listing = dir(fullfile(PCdir,'3d*out'));
   [filehead,data] = read_data(...
         fullfile(listing(1).folder,listing(1).name),...
         'verbose',false);
   npict = numel(listing); 
end

% Maybe write a function?
ne_ = strcmpi('rhos0',filehead.wnames);
ni_ = strcmpi('rhos1',filehead.wnames);
bx_ = strcmpi('bx',filehead.wnames);
by_ = strcmpi('by',filehead.wnames);
bz_ = strcmpi('bz',filehead.wnames);
uex_ = strcmpi('uxs0',filehead.wnames);
uey_ = strcmpi('uys0',filehead.wnames);
uez_ = strcmpi('uzs0',filehead.wnames);
uix_ = strcmpi('uxs1',filehead.wnames);
uiy_ = strcmpi('uys1',filehead.wnames);
uiz_ = strcmpi('uzs1',filehead.wnames);
pe_ = strcmpi('ps0',filehead.wnames);
pi_ = strcmpi('ps1',filehead.wnames);


% Create video
v = VideoWriter(Vname);
v.FrameRate = vFrameRate;
v.open

% create new figure with specified size
hfig = figure(4);
set(hfig,'position', [10, 10, 800, 550]) 
colormap(jet);

for ipict=1:npict
   fprintf('ipict=%d\n',ipict)

   if IsGatheredFile
      [filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);
   else
      [filehead,data] = read_data(...
         fullfile(listing(ipict).folder,listing(ipict).name),...
         'verbose',false);
   end   
   
   data = data.file1;
   x = data.x(:,:,:,1);
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);
   
   ne = data.w(:,:,:,ne_)*1e6;    % [/m^3]
   ni = data.w(:,:,:,ni_)*1e6;    % [/m^3]
   bx   = data.w(:,:,:,bx_);      % [nT]
   by   = data.w(:,:,:,by_);
   bz   = data.w(:,:,:,bz_);
   uex  = data.w(:,:,:,uex_);      % [km/s]
   uey  = data.w(:,:,:,uey_);
   uez  = data.w(:,:,:,uez_);
   uix  = data.w(:,:,:,uix_);
   uiy  = data.w(:,:,:,uiy_);
   uiz  = data.w(:,:,:,uiz_);
   pe   = data.w(:,:,:,pe_);     % [nPa]
   pi   = data.w(:,:,:,pi_);
   
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
   caxis([-0.6 0.6])
   title('U_{iL}')
   
   subplot_tight(4,3,2);
   contourf(yq,zq,uiM./VA,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-0.5 0.5])
   ylabel('z [R_G]');
   title('U_{iM}')
       
   subplot_tight(4,3,3);
   contourf(yq,zq,uiN./VA,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-0.5 0.5])
   title('U_{iN}')
   
   % ueL
   subplot_tight(4,3,4);
   contourf(yq,zq,ueL./VA,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-3 3])
   title('U_{eL}')
   
   % ueM
   subplot_tight(4,3,5);
   contourf(yq,zq,ueM./VA,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0 5])
   title('U_{eM}')
   
   % ueN
   subplot_tight(4,3,6);
   contourf(yq,zq,ueN./VA,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-1 1])
   title('U_{eN}')
   
   % BL
   subplot_tight(4,3,7);
   contourf(yq,zq,bL,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0 180])
   title('B_L [nT]')
   
   % BM (positive pointing in roughly y direction)
   subplot_tight(4,3,8);
   contourf(yq,zq,bM,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-80 80])
   title('B_M [nT]')
   
   % BN (positive pointing upstream)
   subplot_tight(4,3,9);
   contourf(yq,zq,bN,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-50 50])
   title('B_N [nT]')
   
   subplot_tight(4,3,10);
   contourf(yq,zq,piv,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0 7])
   title('P_{i} [nPa]')
   
   subplot_tight(4,3,11);
   contourf(yq,zq,pev,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0 1])
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('P_{e} [nPa]')
   
   subplot_tight(4,3,12);
   contourf(yq,zq,J,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0 0.2])
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
