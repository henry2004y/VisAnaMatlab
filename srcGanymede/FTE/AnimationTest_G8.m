% Script for capturing FTEs and doing statistics.
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
% Hongyang Zhou, hyzhou@umich.edu  11/02/2017


clear; clc
%%

CPCPz3 = load('../ProcessedData/CPCP_G8_600s_z=3.mat','CPCPt');
CPCPz3 = CPCPz3.CPCPt;

CPCPz2 = load('../ProcessedData/CPCP_G8_600s_z=2.mat','CPCPt');
CPCPz2 = CPCPz2.CPCPt;

CPCPb = load('../ProcessedData/CPCP_G8_600s_b=3.mat','CPCPt');
CPCPb = CPCPb.CPCPt;

time = 1:600;
figure;
plot(time,CPCPb,time,CPCPz2,time,CPCPz3);
legend({'b=3','z=2','z=3','new'})

CPCPbfft = fft(CPCPb);
CPCPbfft(1) = [];
n = length(CPCPbfft);
powerb = abs(CPCPbfft(1:floor(n/2))).^2; % power of first half of transform data
maxfreq = 1/2;                   % maximum frequency
freq = (1:n/2)/(n/2)*maxfreq;    % equally spaced frequency grid
periodb = 1./freq;

CPCPz2fft = fft(CPCPz2);
CPCPz2fft(1) = [];
n = length(CPCPz2fft);
powerz2 = abs(CPCPz2fft(1:floor(n/2))).^2; % power of first half of transform data
maxfreq = 1/2;                   % maximum frequency
freq = (1:n/2)/(n/2)*maxfreq;    % equally spaced frequency grid
periodz2 = 1./freq;

CPCPz3fft = fft(CPCPz3);
CPCPz3fft(1) = [];
n = length(CPCPz3fft);
powerz3 = abs(CPCPz3fft(1:floor(n/2))).^2; % power of first half of transform data
maxfreq = 1/2;                   % maximum frequency
freq = (1:n/2)/(n/2)*maxfreq;    % equally spaced frequency grid
periodz3 = 1./freq;


figure;
plot(periodb,powerb,periodz2,powerz2,periodz3,powerz3);
xlim([0 50]); %zoom in on max power
xlabel('Seconds/Cycle')
ylabel('Power')
legend({'b=3','z=2','z=3','new'})

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
%filename = '../../newPIC/run_G8_newPIC_highC_test/GM/box_var_2_n60000_75717.outs'; 
filename = '../../newPIC/box_G8_600s.outs'; 
npict = 600;
v = VideoWriter('G8_test.avi');
v.FrameRate = 10;
v.open

% create new figure with specified size
hfig = figure(3);
set(hfig,'position', [10, 10, 800, 400]) 
colormap(jet);
%fig = gcf;
%ax = fig.CurrentAxes;

for ipict = 1:npict
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
 
   subplot_tight(2,3,1);
   %subplot(2,3,1);
   contourf(yq,zq,b3,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([-60 60])
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('Bnormal [nT]')
   
   subplot_tight(2,3,2);
   %subplot(2,3,2);
   contourf(yq,zq,jv,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   caxis([0.1 0.7])
   xlabel('y [R_G]'); ylabel('z [R_G]');
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
   
   subplot_tight(2,3,3); 
   %subplot(2,3,3);
   contourf(yq,zq,pev,50,'Linestyle','none'); colorbar;
   axis tight equal
   set(gca,'Xdir','reverse')
   %caxis([2 12]);
   caxis([0.2 4]);
   xlabel('y [R_G]'); ylabel('z [R_G]');
   title('pe [nPa]')
   
   hold on;   
   [centroid,boundaries,FTEcountP] = find_FTE(pev,'VarThreshold',...
      threshold_pe,'AreaThreshold',30);
   
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

   %subplot_tight(2,3,[4,5,6]);
   subplot(2,3,[4,5,6]);
   %plot(CPCPt);
   %plot(time,CPCPb,time,CPCPz2,time,CPCPz3);
   plot(time,CPCPnew);
   xline(ipict,'k');  
   %legend({'b','z=2','z=3'})
   
   
   dim = [0.2 0.05 0.05 0.02];
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
