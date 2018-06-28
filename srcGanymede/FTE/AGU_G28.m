% AGU2017 Ganymede G8 & G28 comparison plots
%
% Hongyang Zhou, hyzhou@umich.edu  12/04/2017

clear;clc
%%

%load ../CPCP_G28_b=2.mat

% inew = 1;
% for itime = 1:numel(time)-1
%    if abs(time(itime)-time(itime+1))>0.5
%       timenew(inew) = time(itime);
%       CPCPnew(inew) = CPCPt(itime);
%       inew = inew + 1;
%    end
%    if itime > 1 && abs(time(itime)-time(itime-1))<0.5
%       continue
%    end
% end
% 
% time = timenew;
% CPCPt = CPCPnew;
% 
% clearvars timenew CPCPnew inew

load ../G28_reconnection_rate.mat


%% Find boundary points from steady state solution
filename = '../../run_status_test/3d_idl_G28.out'; % 3d GM outputs
s = 0.5; % compact boundary factor [0,1]

[x3bc,y3bc,z3bc] = find_boundary_points( filename,s );

%% Fit the closed field line boundary

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );


% Plot fit with data.
figure( 'Name', 'poly5' );
h = plot( fitresult );
legend( h, 'poly5, x=x(y,z)', 'Location', 'NorthEast' );
% Label axes
xlabel('x3bc [R_G]')
ylabel('y3bc [R_G]')
zlabel('z3bc [R_G]')
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
ymin = -1.1+5/30; ymax = 1.1-5/30; zmin = -0.5; zmax = 0.75-3/30;
dy = 1/30; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

% Try to rotate the ndgrid by 11 degrees
% Define rotation angle
theta = 11/180*pi; 
% Define rotation matrix
rot = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 


temp = [yq(:),zq(:)]*rot.' ;
sz = numel(yq);
Yrot = reshape(temp(:,1),[sz,1]);
Zrot = reshape(temp(:,2),[sz,1]);

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


%% Get the normal vectors of local coordinates for rotated ndgrid
z = [0 0 1]; % z-direction unit vector
% Initialize local vectors
d1 = Inf(3,size(Xrot,1));
d2 = d1; d3 = d1;

% 
for iP=1:size(Xrot,1)
   d1(:,iP) = cross([Xrot(iP) Yrot(iP) Zrot(iP)],z);
   d2(:,iP) = cross(d1(:,iP),[Xrot(iP) Xrot(iP) Zrot(iP)]);
   d3(:,iP) = [Xrot(iP) Yrot(iP) Zrot(iP)];
   
   % Normalization
   d1(:,iP) = d1(:,iP) / norm(d1(:,iP));
   d2(:,iP) = d2(:,iP) / norm(d2(:,iP));
   d3(:,iP) = d3(:,iP) / norm(d3(:,iP));
end

%% Plots from rotated ndgrid in the local coordinate system

Coef = 1.02; % expansion factor

% Thresholds
threshold_pe = 1.8;
threshold_j = 0.50;
 
filename = '../../newPIC/box_G28_600s.outs'; 
npict = 600; firstpict=221; lastpict=261;
%time = 1:600;
iplot = 1;

% create new figure with specified size
hfig = figure(3);
set(hfig,'position', [10, 10, 800, 600]) 
colormap(jet);


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
   
   bxv= interp3(x, y, z, bx, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:));
   byv= interp3(x, y, z, by, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:)); 
   bzv= interp3(x, y, z, bz, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:));
   jxv= interp3(x, y, z, jx, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:));
   jyv= interp3(x, y, z, jy, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:)); 
   jzv= interp3(x, y, z, jz, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:));  
   pev= interp3(x, y, z, pe, Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:));
   pv = interp3(x, y, z, p,  Coef*Xrot(:), Coef*Yrot(:), Coef*Zrot(:)); 
   
   jv = sqrt(jxv.^2 + jyv.^2 + jzv.^2);
      
   % Transform vectors into local coordinate system
   b1 = Inf(size(Xrot)); b2 = b1; b3 = b1;
      
   % This could potentially be improved!
   for iP=1:size(Xrot,1)
         b1(iP) = [bxv(iP) byv(iP) bzv(iP)]*d1(:,iP);
         b2(iP) = [bxv(iP) byv(iP) bzv(iP)]*d2(:,iP);
         b3(iP) = [bxv(iP) byv(iP) bzv(iP)]*d3(:,iP);
   end
   
   %scatter(Yrot,Zrot,[],b3,'filled'); colorbar; colormap jet; axis tight
   
   b3   = reshape(b3, size(Yq));
   jv   = reshape(jv, size(Yq));
   pev  = reshape(pev,size(Yq));
 
   %subplot_tight(4,3,iplot);
   subaxis(4,3,iplot);
   contourf(Yq,Zq,b3,50,'Linestyle','none');
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
   contourf(Yq,Zq,jv,50,'Linestyle','none');
   
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
   caxis([0.1 0.6])
   
%    hold on;   
%    [centroid,boundaries,FTEcountJ] = find_FTE(jv,'VarThreshold',threshold_j,...
%       'AreaThreshold',40);
%    
%    for k = 1 : FTEcountJ
%       % switch to grid axis
%       thisBoundary = (boundaries{k}-1)*dy + [ymin zmin];
%       thisBoundary = thisBoundary*rot.' ;
%       plot(thisBoundary(:,1), thisBoundary(:,2), 'k', 'LineWidth', 1);
%       % Put the "blob number" labels on the "boundaries" grayscale image.
%       thiscentroid = centroid(k,:)*dy + [zmin ymin];
%       text(thiscentroid(2), thiscentroid(1), num2str(k), ...
%       'FontSize', 10, 'FontWeight', 'Bold')
%    end
%    hold off;
   
   %subplot_tight(4,3,iplot+6); 
   subaxis(4,3,iplot+6);
   contourf(Yq,Zq,pev,50,'Linestyle','none');
   
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
   caxis([0.2 3]);
   
%    hold on;   
%    [centroid,boundaries,FTEcountP] = find_FTE(pev,'VarThreshold',...
%       threshold_pe,'AreaThreshold',30);
%    
%    for k = 1 : numel(boundaries)
%       % switch to grid axis
%       thisBoundary = (boundaries{k}-1)*dy + [ymin zmin];
%       thisBoundary = thisBoundary*rot.' ;
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
%subplot(2,3,[4,5,6]);

plot(time,reconnect_rate_G28,'k','linewidth',1.2);
ylim([0.2 0.4]);
vline(firstpict,'b',sprintf('%ds',firstpict));
vline(firstpict+20,'b',sprintf('%ds',firstpict+20));
vline(lastpict,'b',sprintf('%ds',firstpict+40));
xlabel('time [s]'); 
ylabel('global reconnection rate')
set(gca,'fontsize',14)

% dim = [0.2 0.113 0.05 0.02];
% str = sprintf('it=%d, time=%.1fs, countJ=%d,countP=%d',...
%    filehead.it,filehead.time, FTEcountJ,FTEcountP);
% a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
%    'FontWeight','bold','EdgeColor','none');