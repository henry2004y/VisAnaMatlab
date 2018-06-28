% Script for capturing FTEs and doing statistics in G8 flyby simulation.
% Procedures:
% 1. Get 3D cell-centered status and coordinates from steady state run(t=0)
% 2. Collect all the points with status==3 (which are located on closed 
%    field lines) and find the boundary points at the dayside
% 3. Do interpolation of variable values onto these boundary points for
%    each snapshot. Make contour plots on the projection y-z plane and then
%    animate. Probably you want to pick the max from a shell instead of
%    just a surface cut.
% 4. Capture the edge of peaks and count.
%
% Pay attention to the artificial thresholds I picked. These may be changed
% into an automatic way of picking values.
%
% colormap jet instead of default
%
% Hongyang Zhou, hyzhou@umich.edu  10/12/2017

clear; clc
%% Find boundary points from steady state solution
filename = '~/Ganymede/newPIC/run_G8_newPIC/3d_G8_steady.outs'; % 3d GM outputs
s = 0.5; % compact boundary factor [0,1]

[x3bc,y3bc,z3bc] = find_boundary_points( filename,s );

%% Fit the closed field line boundary with paraboloid

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

% Plot fit with data.
figure( 'Name', 'poly5' );
h = plot( fitresult);
legend( h, 'poly5', 'x3bc vs. y3bc, z3bc', 'Location', 'NorthEast' );
% Label axes
xlabel y3bc
ylabel z3bc
zlabel x3bc
grid on

xx = get(h, 'XData');
yy = get(h, 'YData');
zz = get(h, 'Zdata');
set(h, 'XData', zz, 'YData', xx, 'ZData', yy);

%hold on;
%scatter3(x3bc,y3bc,z3bc,'.'); hold off
%axis tight
xlim([-2 0]); ylim([-1.5 1.5]); zlim([-0.6 0.8]);

%% Generate mesh points from fitted surface
ymin = -1.1+1/15; ymax = 1.1-1/15; zmin = -0.54+1/15; zmax = 0.8-1/15;
dy = 1/30; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

xq = fitresult(yq,zq);


%% Original output variables as a movie
% bx,by,bz,jx,jy,jz,j,b,p

% box outputs
filename = '~/Ganymede/newPIC/run_G8_newPIC/box_FTE_G8_1200s.outs'; 
[~,~,fileinfo] = read_data(filename,'verbose',false);
npict = fileinfo.npictinfiles;

% Parameters and thresholds
Coef = 1.05; % expansion factor
threshold_p = 6.5;
threshold_j = 0.52;

% Create video
v = VideoWriter('../../newPIC/FTE_G8.avi');
v.FrameRate = 10;
v.open

% Create figure with specified size
hfig = figure(3);
set(hfig,'position', [10, 10, 900, 520]);
colormap(jet);

% Loop over snapshots
%1:npict
for ipict = 371:371
   [filehead,data] = read_data(filename,'verbose',false,'npict',ipict);
   
   x = data.file1.x(:,:,:,1);
   y = data.file1.x(:,:,:,2);
   z = data.file1.x(:,:,:,3);
   
   bx = data.file1.w(:,:,:,5);
   by = data.file1.w(:,:,:,6);
   bz = data.file1.w(:,:,:,7);
   b  = sqrt(bx.^2 + by.^2 + bz.^2); 
   p  = data.file1.w(:,:,:,10);
   jx = data.file1.w(:,:,:,11);
   jy = data.file1.w(:,:,:,12);
   jz = data.file1.w(:,:,:,13);
   j  = sqrt(jx.^2 + jy.^2 + jz.^2); 
   %status = data.file1.w(:,:,:,14);
   
   % From ndgrid to meshgrid format
   x  = permute(x,[2 1 3]);
   y  = permute(y,[2 1 3]);
   z  = permute(z,[2 1 3]);
   bx = permute(bx,[2 1 3]);
   by = permute(by,[2 1 3]);
   bz = permute(bz,[2 1 3]);
   b  = permute(b,[2 1 3]);
   p  = permute(p,[2 1 3]);
   jx = permute(jx,[2 1 3]);
   jy = permute(jy,[2 1 3]);
   jz = permute(jz,[2 1 3]);   
   j  = permute(j,[2 1 3]);
   
   bxv= interp3(x, y, z, bx, Coef*xq(:), Coef*yq(:), Coef*zq(:));
   byv= interp3(x, y, z, by, Coef*xq(:), Coef*yq(:), Coef*zq(:)); 
   bzv= interp3(x, y, z, bz, Coef*xq(:), Coef*yq(:), Coef*zq(:)); 
   bv = interp3(x, y, z, b,  Coef*xq(:), Coef*yq(:), Coef*zq(:));
   pv = interp3(x, y, z, p,  Coef*xq(:), Coef*yq(:), Coef*zq(:));
   jxv= interp3(x, y, z, jx, Coef*xq(:), Coef*yq(:), Coef*zq(:));
   jyv= interp3(x, y, z, jy, Coef*xq(:), Coef*yq(:), Coef*zq(:)); 
   jzv= interp3(x, y, z, jz, Coef*xq(:), Coef*yq(:), Coef*zq(:));   
   jv = interp3(x, y, z, j,  Coef*xq(:), Coef*yq(:), Coef*zq(:));
 
   % Note that there are still some NaN on the boundaries! I could, for
   % example, reduce the meshgrid size by 4 in each dimension.
   % Now there are none.
   
   bxv = reshape(bxv,size(xq));
   byv = reshape(byv,size(xq));
   bzv = reshape(bzv,size(xq));
   jxv = reshape(jxv,size(xq));
   jyv = reshape(jyv,size(xq));
   jzv = reshape(jzv,size(xq));
   bv  = reshape(bv,size(xq));
   jv  = reshape(jv,size(xq));
   pv  = reshape(pv,size(xq));
      
   subplot_tight(3,3,1);
   %scatter(y3bc,z3bc,[],bxv,'filled'); colorbar; colormap jet; axis tight
   
   contourf(yq,zq,bxv,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-100 100]);
   title('Bx [nT]');

   subplot_tight(3,3,2);
   %scatter(y3bc,z3bc,[],byv,'filled'); colorbar; colormap jet; axis tight
   contourf(yq,zq,byv,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')   
   caxis([-60 60]);
   title('By [nT]');
   
   subplot_tight(3,3,3);
   %scatter(y3bc,z3bc,[],bzv,'filled'); colorbar;
   contourf(yq,zq,bzv,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')   
   caxis([-80 80]);
   title('Bz [nT]');
   
   subplot_tight(3,3,4);
   %scatter(y3bc,z3bc,[],jxv,'filled'); colorbar;
   contourf(yq,zq,jxv,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')   
   caxis([-0.3 0.2]);
   title('jx');
   
   subplot_tight(3,3,5);
   %scatter(y3bc,z3bc,[],jyv,'filled'); colorbar;
   contourf(yq,zq,jyv,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-0.8 0.0]);
   title('jy');
   
   subplot_tight(3,3,6);
   %scatter(y3bc,z3bc,[],jzv,'filled'); colorbar;
   contourf(yq,zq,jzv,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')   
   caxis([-0.4 0.4]);
   title('jz');
   
   subplot_tight(3,3,7);
   %scatter(y3bc,z3bc,[],jv,'filled'); colorbar;
   contourf(yq,zq,jv,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')   
   caxis([0.1 0.9])
   title('j')
   
   hold on;   
   [centroid,boundaries,FTEcountJ] = find_FTE(jv,'VarThreshold',threshold_j,...
      'AreaThreshold',40);
   
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
   
   subplot_tight(3,3,8);
   %scatter(y3bc,z3bc,[],bv,'filled'); colorbar;
   contourf(yq,zq,bv,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')   
   caxis([20 140]);
   title('B [nT]')
   
   subplot_tight(3,3,9);
   %scatter(y3bc,z3bc,[],pv,'filled'); colorbar;
   contourf(yq,zq,pv,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')   
   caxis([2 12]);
   title('p [nPa]')
   
   hold on;   
   [centroid,boundaries,FTEcountP] = find_FTE(pv,'VarThreshold',threshold_p,...
      'AreaThreshold',40);
   
   %pthres(ipict) = graythresh((pv-min(pv(:))) ./ (max(pv(:))-min(pv(:))));
   
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
   
   
   dim = [0.1 0.01 0.1 0.045];
   str = sprintf('it=%d, time=%.1fs, countJ=%d,countP=%d',...
      filehead.it,filehead.time, FTEcountJ,FTEcountP);
   a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
      'FontWeight','bold','EdgeColor','none');
   
   frame = getframe(gcf);
   
   set(gca,'nextplot','replacechildren');
   
   clf;
   writeVideo(v,frame);
   
end

v.close
close(3)
