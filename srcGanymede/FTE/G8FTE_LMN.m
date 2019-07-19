% Script for capturing and analyzing FTEs in G8 flyby simulation.
% Procedures:
% 1. Get 3D cell-centered status and coordinates from steady state run(t=0)
% 2. Collect all the points with status==3 (which are located on closed 
%    field lines) and find the boundary points at the upstream.
% 3. Interpolate variables onto these boundary points for each snapshot.
%    Make contour plots in the boundary-normal coordinates and animate. The
%    boundary-normal coordinates is defined as follows: one axis (d3) 
%    normal to the boundary surface pointing outward(upstream); one axis 
%    (d1) tangent to the surface as \vec{d1}=\vec{d3}x\vec{z}; the other 
%    (d2) complete the right-hand rule. Note that in the following code, 
%    I only care about the orientation, so d1, d2 and d3 may not strictly
%    follow the right-hand rule. However, it is not needed to be truly
%    right-handed as long as you know the direction of vectors.
%    Probably you want to pick the max value from a shell instead of single
%    surface cut. (Gabor suggested it, but I`ve never tried.) 
% 4. Capture the edge of peaks and count. I tried some technics from image
%    processing, but the results are not good enough. The number of FTEs is
%    hard to count, because we don`t have a valid definition of FTE.
%
% Pay attention to the artificial thresholds I picked. These may be changed
% into an automatic way of picking values.
%
% In the coutour plots, the y axis is reversed because I want to show a
% view from upstream to downstream in GPhiO coordinate system.
%
% In each processing step, there are some coefficients that can be
% modified. Be careful for those numbers!
%
% Hongyang Zhou, hyzhou@umich.edu  11/02/2017, version 1.1
%
% modified 01/03/2018, version 1.2
% modified 06/29/2018, version 1.3

clear; clc
%% Parameters
% 3D GM outputs
filename = '~/Documents/research/Ganymede/SteadyRun/3d_G8_idl.out';
s = 0.5; % compact boundary factor [0,1]
DoPlot = true;
xThres = 1.5;
rThres = -1.1;

% GM upstream box outputs
fnameFTE = '~/Ganymede/newPIC/run_G8_newPIC/box_FTE_G8_1200s.outs';
% Output movie
Vname = 'test.avi';
vFrameRate = 10;

% Criteria for surface contour and FTE identification
Coef         = 1.02; % expansion factor from original surface fit 
threshold_pe = 2.1;
threshold_j  = 0.52; 

%% Find boundary points from steady state solution

[x3bc,y3bc,z3bc] = find_boundary_points( filename,s,DoPlot,xThres,rThres );

%% Fit the closed field line boundary with hypersurface fit

[fitresult,gof] = surface_fit(x3bc,y3bc,z3bc);

%% Generate mesh points from fitted surface
ymin = -1.1+1/15; ymax = 1.1-1/15; zmin = -0.54+1/15; zmax = 0.8-1/15;
dy = 1/30; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

xq = fitresult(yq,zq);

% Transform into local coordinate system

% Calculate the normal direction to the fitted surface
[V, W] = differentiate(fitresult, yq, zq);

U = -ones(size(V));

% [U,V,W]
figure(1); hold on
quiver3(xq,yq,zq,U,V,W,2,'color','r')
hold off

% get the three local directions
% dipole-direction unit vector
unitDipole = [18 -51.82 716.8]/sqrt(18^2+51.82^2+716.8^2);
% Initialize local vectors
dL = Inf(3,size(xq,1),size(xq,2));
dM = dL; dN = dL;

% This part could potentially be optimized!
for ix=1:size(xq,1); for iy=1:size(xq,2)
   dN(:,ix,iy) = [U(ix,iy) V(ix,iy) W(ix,iy)];
   dM(:,ix,iy) = cross(dN(:,ix,iy),unitDipole);
   dL(:,ix,iy) = cross(dM(:,ix,iy),dN(:,ix,iy));
      
   % Normalization
   dL(:,ix,iy) = dL(:,ix,iy) / norm(dL(:,ix,iy));
   dM(:,ix,iy) = dM(:,ix,iy) / norm(dM(:,ix,iy));
   dN(:,ix,iy) = dN(:,ix,iy) / norm(dN(:,ix,iy));
end; end

% Test this for normalization
%dN = Inf(3,size(xq,1)*size(xq,2));
% dN = [U(:) V(:) W(:)];
% dN = bsxfun(@rdivide,dN,sum(dN,2));
% dN = reshape(dN,size(xq));
%bsxfun(@times,)

%% Visualizing in local coordinates as a movie

[~,~,fileinfo] = read_data(fnameFTE,'verbose',false);
npict = fileinfo.npictinfiles;

FTEcount = zeros(npict,2); % # counts for FTE with J and P respectively

% Create video
v = VideoWriter(Vname);
v.FrameRate = vFrameRate;
v.open

% Create figure with specified size
hfig = figure(3);
set(hfig,'position', [10, 10, 800, 300]) 
colormap(jet);

% Loop over snapshots
for ipict=1:1%npict
   fprintf('ipict=%d\n',ipict)
   [filehead,data] = read_data(fnameFTE,'verbose',false,'npict',ipict);
   
   data = data.file1;
   
   x = data.x(:,:,:,1);
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);

   ux = data.w(:,:,:,2);
   uy = data.w(:,:,:,3);
   uz = data.w(:,:,:,4);   
   bx = data.w(:,:,:,5);
   by = data.w(:,:,:,6);
   bz = data.w(:,:,:,7);
   pe = data.w(:,:,:,9);
   p  = data.w(:,:,:,10);
   jx = data.w(:,:,:,11);
   jy = data.w(:,:,:,12);
   jz = data.w(:,:,:,13);
   %j  = sqrt(jx.^2 + jy.^2 + jz.^2); 
   
   % From ndgrid to meshgrid format
   x  = permute(x,[2 1 3]);
   y  = permute(y,[2 1 3]);
   z  = permute(z,[2 1 3]);
   ux = permute(ux,[2 1 3]);
   uy = permute(uy,[2 1 3]);
   uz = permute(uz,[2 1 3]);   
   bx = permute(bx,[2 1 3]);
   by = permute(by,[2 1 3]);
   bz = permute(bz,[2 1 3]);
   jx = permute(jx,[2 1 3]);
   jy = permute(jy,[2 1 3]);
   jz = permute(jz,[2 1 3]); 
   pe = permute(pe,[2 1 3]);
   p  = permute(p,[2 1 3]);
   %j  = permute(j,[2 1 3]);

   uxv= interp3(x, y, z, ux, Coef*xq, Coef*yq, Coef*zq);
   uyv= interp3(x, y, z, uy, Coef*xq, Coef*yq, Coef*zq); 
   uzv= interp3(x, y, z, uz, Coef*xq, Coef*yq, Coef*zq);   
   bxv= interp3(x, y, z, bx, Coef*xq, Coef*yq, Coef*zq);
   byv= interp3(x, y, z, by, Coef*xq, Coef*yq, Coef*zq); 
   bzv= interp3(x, y, z, bz, Coef*xq, Coef*yq, Coef*zq);
   jxv= interp3(x, y, z, jx, Coef*xq, Coef*yq, Coef*zq);
   jyv= interp3(x, y, z, jy, Coef*xq, Coef*yq, Coef*zq); 
   jzv= interp3(x, y, z, jz, Coef*xq, Coef*yq, Coef*zq);  
   pev= interp3(x, y, z, pe, Coef*xq, Coef*yq, Coef*zq);
   pv = interp3(x, y, z, p,  Coef*xq, Coef*yq, Coef*zq); 
   %jv = interp3(x, y, z, j,  Coef*xq, Coef*yq, Coef*zq);
 
   % Note that there may be some NaN because they are out of bound.
   % You could reduce the meshgrid size if that happens.
   
   jv = sqrt(jxv.^2 + jyv.^2 + jzv.^2);
      
   % Transform vectors into local coordinate system
   uL = Inf(size(xq)); uM = uL; uN = uL;
   bL = Inf(size(xq)); bM = bL; bN = bL;
   jL = Inf(size(xq)); jM = jL; jN = jL; 
      
   % This could potentially be improved!
   for ix=1:size(xq,1)
      for iy=1:size(xq,2)
         bL(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
         bM(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
         bN(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dN(:,ix,iy);
         jL(ix,iy) = [jxv(ix,iy) jyv(ix,iy) jzv(ix,iy)]*dL(:,ix,iy);
         jM(ix,iy) = [jxv(ix,iy) jyv(ix,iy) jzv(ix,iy)]*dM(:,ix,iy);
         jN(ix,iy) = [jxv(ix,iy) jyv(ix,iy) jzv(ix,iy)]*dN(:,ix,iy);         
      end
   end
 
   % Bn (positive pointing upstream)
   subplot_tight(2,3,1);   
   contourf(yq,zq,bN,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-60 60])
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('B_N [nT]')

   % Bm (positive pointing in roughly y direction)
   subplot_tight(2,3,2);   
   contourf(yq,zq,bM,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-60 30])
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('B_M [nT]')   
   
   % Bl
   subplot_tight(2,3,3);
   contourf(yq,zq,bL,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   %caxis([0.1 0.7])
   
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('B_L')
   
   hold on;   

   % Choose the current threshold based on distribution
   [mean,dev] = normfit(jv(:));
   if dev/mean>=0.3
      [centroid,boundaries,FTEcountJ] = find_FTE(jv,'VarThreshold',...
         mean+1.5*dev,'AreaThreshold',40);
   else
      [centroid,boundaries,FTEcountJ] = find_FTE(jv,'VarThreshold',...
         threshold_j,'AreaThreshold',40);      
   end

   % Plot the identified FTE boundary onto the plane
   for k = 1 : FTEcountJ
      % switch to grid axis
      thisBoundary = (boundaries{k}-1)*dy + [ymin zmin];
      plot(thisBoundary(:,1), thisBoundary(:,2), 'k', 'LineWidth', 1);
      % Put the "blob number" labels on the "boundaries" grayscale image.
      thiscentroid = centroid(k,:)*dy + [zmin ymin];
      text(thiscentroid(2), thiscentroid(1), num2str(k), ...
      'FontSize', 10, 'FontWeight', 'Bold')
   end
   hold off;
   
   % jL (positive in rougly z direction)
   subplot_tight(2,3,4);   
   contourf(yq,zq,jL,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   %caxis([-200 200])
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('j_L [\mu A/m^2]')

   % jM (positive pointing in roughly y direction)
   subplot_tight(2,3,5);   
   contourf(yq,zq,jM,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   %caxis([-100 100])
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('j_M [\mu A/m^2]') 
   
   
   % jN
   subplot_tight(2,3,6); 
   contourf(yq,zq,jN,'Linestyle','none'); colorbar; 
   axis tight equal
   set(gca,'Xdir','reverse')
   %caxis([0.2 4]);
   %caxis([2 12]);
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('j_N')
   
   hold on;   
 
   % Choose the current threshold based on distribution   
   [mean,dev] = normfit(pev(:));
   if dev/mean>=0.38
      [centroid,boundaries,FTEcountP] = find_FTE(pev,'VarThreshold',...
         mean+1.5*dev,'AreaThreshold',30);
   else
      [centroid,boundaries,FTEcountP] = find_FTE(pev,'VarThreshold',...
         threshold_pe,'AreaThreshold',30);      
   end

   % Plot the identified FTE boundary onto the plane   
   for k = 1 : numel(boundaries)
      % switch to grid axis
      thisBoundary = (boundaries{k}-1)*dy + [ymin zmin];
      plot(thisBoundary(:,1), thisBoundary(:,2), 'k', 'LineWidth', 1);
      % Put the "blob number" labels on the "boundaries" grayscale image.
      thiscentroid = centroid(k,:)*dy + [zmin ymin];
      text(thiscentroid(2), thiscentroid(1), num2str(k), ...
      'FontSize', 10, 'FontWeight', 'Bold')
   end
   hold off;
   
   dim = [0.17 0.05 0.05 0.02];
   str = sprintf('it=%d, time=%.1fs, countJ=%d,countP=%d',...
      filehead.it,filehead.time, FTEcountJ,FTEcountP);
   a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
      'FontWeight','bold','EdgeColor','none');
   
   frame = getframe(gcf);
   
   % remove the plots while keeping axis properties
   set(gca,'nextplot','replacechildren');
   
   % clear current figure window and write to video
   clf;
   writeVideo(v,frame);
   
   % Record the counts identified by image processing algorithm
   FTEcount(ipict,1) = FTEcountJ;
   FTEcount(ipict,2) = FTEcountP;

end

v.close
close(3)
