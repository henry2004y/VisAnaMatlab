% AGU2017 Ganymede G8 & G28 comparison plots
%
% Hongyang Zhou, hyzhou@umich.edu  12/04/2017

clear;clc
%%

%load ../CPCP_G8_new_x=1.mat
%load ../CPCP_G8_600s_z=2.mat
load ../G8_reconnection_rate.mat

%% Find boundary points from steady state solution
filename = '../../run_status_test/GM/3d_idl.outs'; % 3d GM outputs
s = 0.5; % compact boundary factor [0,1]

[x3bc,y3bc,z3bc] = find_boundary_points( filename,s );

%% Fit the closed field line boundary

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

%% Generate mesh points from fitted surface
ymin = -1.1+1/15; ymax = 1.1-1/15; zmin = -0.54+1/15; zmax = 0.8-1/15;
dy = 1/30; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

xq = fitresult(yq,zq);

%% get the three local directions
z = [0 0 1]; % z-direction unit vector
% Initialize local vectors
d1 = Inf(3,size(xq,1),size(xq,2));
d2 = d1; d3 = d1;

% This part could potentially be optimized!
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      d1(:,ix,iy) = cross([xq(ix,iy) yq(ix,iy) zq(ix,iy)],z);
      d2(:,ix,iy) = cross(d1(:,ix,iy),[xq(ix,iy) yq(ix,iy) zq(ix,iy)]);
      d3(:,ix,iy) = [xq(ix,iy) yq(ix,iy) zq(ix,iy)];
      
      % Normalization
      d1(:,ix,iy) = d1(:,ix,iy) / norm(d1(:,ix,iy));
      d2(:,ix,iy) = d2(:,ix,iy) / norm(d2(:,ix,iy));
      d3(:,ix,iy) = d3(:,ix,iy) / norm(d3(:,ix,iy));
   end
end

%% Output variables in the local coordinate system as a movie

Coef = 1.03; % expansion factor

% Thresholds
threshold_pe = 2.0;
threshold_j = 0.52; 
filename = '../../newPIC/box_G8_600s.outs'; 
npict = 600; firstpict = 110; lastpict = 150;
time = 1:600;

% create new figure with specified size
hfig = figure;
set(hfig,'position', [10, 10, 800, 600]) 
colormap(jet);
%fig = gcf;
%ax = fig.CurrentAxes;
iplot = 1;

%[ha, pos] = tight_subplot(3,3,[.01 .03],[.1 .01],[.01 .01]);

for ipict = firstpict:20:lastpict
   [filehead,data] = read_data(filename,'verbose',false,'npict',ipict);
   
   x = data.file1.x(:,:,:,1);
   y = data.file1.x(:,:,:,2);
   z = data.file1.x(:,:,:,3);
   
   bx = data.file1.w(:,:,:,5);
   by = data.file1.w(:,:,:,6);
   bz = data.file1.w(:,:,:,7);
   pe = data.file1.w(:,:,:,9);
   p  = data.file1.w(:,:,:,10);
   jx = data.file1.w(:,:,:,11);
   jy = data.file1.w(:,:,:,12);
   jz = data.file1.w(:,:,:,13);
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
   
   bxv= interp3(x, y, z, bx, Coef*xq(:), Coef*yq(:), Coef*zq(:));
   byv= interp3(x, y, z, by, Coef*xq(:), Coef*yq(:), Coef*zq(:)); 
   bzv= interp3(x, y, z, bz, Coef*xq(:), Coef*yq(:), Coef*zq(:));
   jxv= interp3(x, y, z, jx, Coef*xq(:), Coef*yq(:), Coef*zq(:));
   jyv= interp3(x, y, z, jy, Coef*xq(:), Coef*yq(:), Coef*zq(:)); 
   jzv= interp3(x, y, z, jz, Coef*xq(:), Coef*yq(:), Coef*zq(:));  
   pev= interp3(x, y, z, pe, Coef*xq(:), Coef*yq(:), Coef*zq(:));
   pv = interp3(x, y, z, p,  Coef*xq(:), Coef*yq(:), Coef*zq(:)); 
   %jv = interp3(x, y, z, j,  Coef*xq(:), Coef*yq(:), Coef*zq(:));
 
   % Note that there are still some NaN on the boundaries! I could, for
   % example, reduce the meshgrid size by 4 in each dimension.
   % Now there are none.
   
   bxv = reshape(bxv,size(xq));
   byv = reshape(byv,size(xq));
   bzv = reshape(bzv,size(xq));
   jxv = reshape(jxv,size(xq));
   jyv = reshape(jyv,size(xq));
   jzv = reshape(jzv,size(xq));
   %jv  = reshape(jv,size(xq));
   pev = reshape(pev,size(xq));
   pv  = reshape(pv,size(xq));
   
   jv = sqrt(jxv.^2 + jyv.^2 + jzv.^2);
      
   % Transform vectors into local coordinate system
   u1 = Inf(size(xq)); u2 = u1; u3 = u1;
   b1 = Inf(size(xq)); b2 = b1; b3 = b1;
      
   % This could potentially be improved!
   for ix=1:size(xq,1)
      for iy=1:size(xq,2)
         b1(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*d1(:,ix,iy);
         b2(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*d2(:,ix,iy);
         b3(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*d3(:,ix,iy);
      end
   end
 
   %subplot_tight(4,3,iplot);
   subaxis(4,3,iplot);
   %subplot(4,3,iplot);
   %axes(ha(iplot));
   
   contourf(yq,zq,b3,50,'Linestyle','none'); 
   title(sprintf('t=%ds',ipict));
   
   if iplot == 1
      text(1.2,0,'Bnormal [nT]','rotation',90,...
         'fontsize',14,'FontWeight','bold',...
         'horizontalalignment','center','verticalalignment','bottom');
   elseif iplot == 3
      % Get the current axis size
      originalSize = get(gca, 'Position');
      colorbar; 
      set(gca,'Position', originalSize);
   end
   
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-60 60])
   
   %subplot_tight(4,3,iplot+3);
   subaxis(4,3,iplot+3);
   %subplot(4,3,iplot+3);
   %axes(ha(iplot+3));
   
   contourf(yq,zq,jv,50,'Linestyle','none'); 
   
   if iplot == 1
      text(1.2,0,'j [\mu A/m^2]','rotation',90,...
         'fontsize',14,'FontWeight','bold',...
         'horizontalalignment','center','verticalalignment','bottom');
   elseif iplot == 3
      % Get the current axis size
      originalSize = get(gca, 'Position');
      c = colorbar;
      set(gca,'Position', originalSize);
   end
   
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0.1 0.7])
   
%    hold on;   
%    [centroid,boundaries,FTEcountJ] = find_FTE(jv,'VarThreshold',threshold_j,...
%       'AreaThreshold',40);
%    
%    for k = 1 : FTEcountJ
%       % switch to grid axis
%       thisBoundary = (boundaries{k}-1)*dy + [ymin zmin];
%       plot(thisBoundary(:,1), thisBoundary(:,2), 'k', 'LineWidth', 1);
%       % Put the "blob number" labels on the "boundaries" grayscale image.
%       thiscentroid = centroid(k,:)*dy + [zmin ymin];
%       text(thiscentroid(2), thiscentroid(1), num2str(k), ...
%       'FontSize', 10, 'FontWeight', 'Bold')
%    end
%    hold off;
   
   %subplot_tight(4,3,iplot+6); 
   subaxis(4,3,iplot+6);
   %subplot(4,3,iplot+6);
   %axes(ha(iplot+6));
   
   contourf(yq,zq,pev,50,'Linestyle','none'); 

   if iplot == 1
      text(1.2,0,'P_e [nPa]','rotation',90,...
         'fontsize',14,'FontWeight','bold',...
         'horizontalalignment','center','verticalalignment','bottom');
   elseif iplot == 2
      xlabel('Y [R_G]','FontWeight','bold','fontsize',14,...
         'horizontalalignment','center','verticalalignment','bottom');
      ylabel('Z [R_G]','FontWeight','bold','fontsize',14,...
         'horizontalalignment','center','verticalalignment','bottom');
   elseif iplot == 3
      % Get the current axis size
      originalSize = get(gca, 'Position');
      colorbar; 
      set(gca,'Position', originalSize);      
   end   
       
   axis tight equal
   set(gca,'Xdir','reverse')
   %caxis([2 12]);
   caxis([0.2 4]);

%    hold on;   
%    [centroid,boundaries,FTEcountP] = find_FTE(pev,'VarThreshold',...
%       threshold_pe,'AreaThreshold',30);
%    
%    for k = 1 : numel(boundaries)
%       % switch to grid axis
%       thisBoundary = (boundaries{k}-1)*dy + [ymin zmin];
%       plot(thisBoundary(:,1), thisBoundary(:,2), 'k', 'LineWidth', 1);
%       % Put the "blob number" labels on the "boundaries" grayscale image.
%       thiscentroid = centroid(k,:)*dy + [zmin ymin];
%       text(thiscentroid(2), thiscentroid(1), num2str(k), ...
%       'FontSize', 10, 'FontWeight', 'Bold')
%    end
%    hold off;   
     
   iplot = iplot + 1;
end

%subplot_tight(4,3,[10,11,12]);
subaxis(4,3,[10,11,12]);
%subplot(4,3,[10,11,12]);
%figure

plot(time(1:300),reconnect_rate_G8(1:300),'k','linewidth',1.2);
vline(firstpict,'b',sprintf('%ds',firstpict));
vline(firstpict+20,'b',sprintf('%ds',firstpict+20));
vline(lastpict,'b',sprintf('%ds',firstpict+40));
xlabel('time [s]'); 
ylabel('global reconnection rate')
set(gca,'fontsize',14)


% dim = [0.2 0.05 0.05 0.02];
% str = sprintf('it=%d, time=%.1fs, countJ=%d,countP=%d',...
%    filehead.it,filehead.time, FTEcountJ,FTEcountP);
% a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
%    'FontWeight','bold','EdgeColor','none');