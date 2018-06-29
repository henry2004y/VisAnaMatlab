% Script for capturing and analyzing FTEs in G28 flyby simulation.
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
%    the right-hand rule. However, it is not needed to be truly
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
% Comparing with G8, I rotate the orientation of meshgrid to fit the tilted
% surface due to strong By in the upstream condition.
%
% Hongyang Zhou, hyzhou@umich.edu  11/02/2017, version 1.1
%
% modified 01/03/2018, version 1.2
% modified 06/29/2018, version 1.3


clear; clc
%% Parameters
% 3D GM outputs
filename = '~/Ganymede/newPIC/run_G28_newPIC/3d_G28_steady.outs';
s = 0.5; % compact boundary factor [0,1]
DoPlot = false;
xThres = 1.5;
rThres = -1.1;

% GM upstream box outputs
fnameFTE = '~/Ganymede/newPIC/run_G28_newPIC/box_FTE_G28_1200s.outs';
% Output movie
Vname = 'test.avi';
vFrameRate = 10;

% Criteria for surface contour and FTE identification
Coef         = 1.02; % expansion factor from original surface fit
threshold_pe = 1.8;
threshold_j  = 0.50;

%% Find boundary points from steady state solution

[x3bc,y3bc,z3bc] = find_boundary_points( filename,s,DoPlot,xThres,rThres );

%% Fit the closed field line boundary with hypersurface

[fitresult,gof] = surface_fit(x3bc,y3bc,z3bc);

%% Generate mesh points from fitted surface
ymin = -1.1+5/30; ymax = 1.1-5/30; zmin = -0.5; zmax = 0.75-3/30;
dy = 1/30; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

% Try to rotate the ndgrid by 13 degrees
% Define rotation angle
theta = 13/180*pi; 
% Define rotation matrix
rot = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 

temp = [yq(:),zq(:)]*rot.' ;
Yrot = reshape(temp(:,1),[numel(yq),1]);
Zrot = reshape(temp(:,2),[numel(yq),1]);

Xrot = fitresult(Yrot,Zrot);

Xq = reshape(Xrot,size(yq));
Yq = reshape(Yrot,size(yq));
Zq = reshape(Zrot,size(yq));

%% Transform into local coordinate system

% Try to do the same for rotated ndgrid points
[V, W] = differentiate(fitresult, Yrot, Zrot);

U = -ones(size(V));

figure(2); hold on
quiver3(Xrot,Yrot,Zrot,U,V,W,2,'color','r')
xlabel('x'); ylabel('y'); zlabel('z');
hold off; box on

% Get the normal vectors of local coordinates for rotated ndgrid
unitz = [0 0 1]; % z-direction unit vector
% Initialize local vectors
dL = Inf(3,size(Xrot,1));
dM = dL; dN = dL;

for iP=1:size(Xrot,1)
   dN(:,iP) = [U(iP) V(iP) W(iP)];
   dL(:,iP) = cross(dN(:,iP),unitz);
   dM(:,iP) = cross(dL(:,iP),dN(:,iP));

   % Normalization
   dL(:,iP) = dL(:,iP) / norm(dL(:,iP));
   dM(:,iP) = dM(:,iP) / norm(dM(:,iP));
   dN(:,iP) = dN(:,iP) / norm(dN(:,iP));
end

%% Visualizing in local coordinates as a movie

[~,~,fileinfo] = read_data(fnameFTE,'verbose',false);
npict = fileinfo.npictinfiles;

FTEcount = zeros(npict,2); % # counts for FTE with J and P respectively

% Create video
v = VideoWriter(Vname);
v.FrameRate = vFrameRate;
v.open

% create new figure with specified size
hfig = figure(3);
set(hfig,'position', [10, 10, 800, 200]) 
colormap(jet);

% Loop over snapshots
for ipict = 1:npict
   [filehead,data] = read_data(fnameFTE,'verbose',false,'npict',ipict);
   
   data = data.file1;
   
   x = data.x(:,:,:,1);
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);
   
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
   bx = permute(bx,[2 1 3]);
   by = permute(by,[2 1 3]);
   bz = permute(bz,[2 1 3]);
   jx = permute(jx,[2 1 3]);
   jy = permute(jy,[2 1 3]);
   jz = permute(jz,[2 1 3]); 
   pe = permute(pe,[2 1 3]);
   p  = permute(p,[2 1 3]);
   %j  = permute(j,[2 1 3]);
   
   bxv= interp3(x, y, z, bx, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:));
   byv= interp3(x, y, z, by, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:)); 
   bzv= interp3(x, y, z, bz, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:));
   jxv= interp3(x, y, z, jx, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:));
   jyv= interp3(x, y, z, jy, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:)); 
   jzv= interp3(x, y, z, jz, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:));  
   pev= interp3(x, y, z, pe, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:));
   pv = interp3(x, y, z, p,  Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:)); 
   %jv = interp3(x, y, z, j,  Coef*xq(:), Coef*yq(:), Coef*zq(:));
 
   % Note that there may be some NaN because they are out of bound.
   % You could reduce the meshgrid size if that happens.
   
%    bxv = reshape(bxv,size(xq));
%    byv = reshape(byv,size(xq));
%    bzv = reshape(bzv,size(xq));
%    jxv = reshape(jxv,size(xq));
%    jyv = reshape(jyv,size(xq));
%    jzv = reshape(jzv,size(xq));
%    %jv  = reshape(jv,size(xq));
%    pev = reshape(pev,size(xq));
%    pv  = reshape(pv,size(xq));
   
   jv = sqrt(jxv.^2 + jyv.^2 + jzv.^2);
      
   % Transform vectors into local coordinate system
   b1 = Inf(size(Xrot)); b2 = b1; b3 = b1;
      
   % This could potentially be improved!
   for iP=1:size(Xrot,1)
         b1(iP) = [bxv(iP) byv(iP) bzv(iP)]*dL(:,iP);
         b2(iP) = [bxv(iP) byv(iP) bzv(iP)]*dM(:,iP);
         b3(iP) = [bxv(iP) byv(iP) bzv(iP)]*dN(:,iP);
   end
   
   %scatter(Yrot,Zrot,[],b3,'filled'); colorbar; colormap jet; axis tight
   
   b3   = reshape(b3, size(Yq));
   jv   = reshape(jv, size(Yq));
   pev  = reshape(pev,size(Yq));
 
   % B normal (positive pointing upstream)   
   subplot(1,3,1);
   contourf(Yq,Zq,b3,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-60 60])
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('Bnormal [nT]')
     
   % Total current density  
   subplot(1,3,2);
   contourf(Yq,Zq,jv,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0.1 0.6])
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('j')
   
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
      thisBoundary = thisBoundary*rot.' ;
      plot(thisBoundary(:,1), thisBoundary(:,2), 'k', 'LineWidth', 1);
      % Put the "blob number" labels on the "boundaries" grayscale image.
      thiscentroid = centroid(k,:)*dy + [zmin ymin];
      text(thiscentroid(2), thiscentroid(1), num2str(k), ...
      'FontSize', 10, 'FontWeight', 'Bold')
   end
   hold off;
   
   % Electron pressure
   subplot(1,3,3); 
   contourf(Yq,Zq,pev,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   %caxis([2 12]);
   caxis([0.2 3]);
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('pe [nPa]')
   
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
      thisBoundary = thisBoundary*rot.' ;
      plot(thisBoundary(:,1), thisBoundary(:,2), 'k', 'LineWidth', 1);
      % Put the "blob number" labels on the "boundaries" grayscale image.
      thiscentroid = centroid(k,:)*dy + [zmin ymin];
      text(thiscentroid(2), thiscentroid(1), num2str(k), ...
      'FontSize', 10, 'FontWeight', 'Bold')
   end
   hold off;
   
   dim = [0.2 0.113 0.05 0.02];
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
