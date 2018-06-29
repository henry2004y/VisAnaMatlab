% Contour plot on a curve surface --- slice
%
% Test on magnetospause surface plot. Later use for FTE_PC_G8
%
% Hongyang Zhou, hyzhou@umich.edu 02/28/2018


clear; clc
%% Find boundary points from steady state solution
filename = '~/Ganymede/newPIC/run_G8_newPIC/3d_G8_steady.outs'; % 3d GM outputs
s = 0.5; % compact boundary factor [0,1]

[x3bc,y3bc,z3bc] = find_boundary_points( filename,s );


%% Fit the closed field line boundary with hypersurface

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

%% Generate mesh points from fitted surface and calculate LMN directions
ymin = -1.5; ymax = 1.5; zmin = -1; zmax = 1;
dy = 1/30; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

xq = fitresult(yq,zq);

% Boundary normal coordinate
unitz = [0 0 1]; % z-direction unit vector
% Initialize local vectors
dM = Inf(3,size(xq,1),size(xq,2)); dL = dM; dN = dM;

% Calculate the normal direction to the fitted surface
[V, W] = differentiate(fitresult, yq, zq);
U = -ones(size(V));


% Get the three local directions
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      dN(:,ix,iy) = [U(ix,iy) V(ix,iy) W(ix,iy)];
      dM(:,ix,iy) = cross(dN(:,ix,iy),unitz);
      dL(:,ix,iy) = cross(dM(:,ix,iy),dN(:,ix,iy));
      
      % Normalization
      dM(:,ix,iy) = dM(:,ix,iy) / norm(dM(:,ix,iy));
      dL(:,ix,iy) = dL(:,ix,iy) / norm(dL(:,ix,iy));
      dN(:,ix,iy) = dN(:,ix,iy) / norm(dN(:,ix,iy));
   end
end



%%
filenamePC='~/Ganymede/newPIC/run_G8_newPIC/3d_fluid_35.outs';
ipict = 10;
[filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);

data = data.file1;
x = data.x(:,:,:,1);       % [Rg]
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);
Rhoe = data.w(:,:,:,1);    % [amu/cc]
Rhoi = data.w(:,:,:,2);    % [amu/cc]
Bx = data.w(:,:,:,3);      % [nT]
By = data.w(:,:,:,4);
Bz = data.w(:,:,:,5);
Ex = data.w(:,:,:,6)*1e-3; % [mV/m]
Ey = data.w(:,:,:,7)*1e-3; % [mV/m]
Ez = data.w(:,:,:,8)*1e-3; % [mV/m]
Uex = data.w(:,:,:,9);     % [km/s]
Uey = data.w(:,:,:,10);
Uez = data.w(:,:,:,11);
Uix = data.w(:,:,:,12);
Uiy = data.w(:,:,:,13);
Uiz = data.w(:,:,:,14);

% The original data is saved in ndgrid format. For streamline and
% isonormals functions, the input should be in meshgrid format.
x  = permute(x,[2 1 3]);
y  = permute(y,[2 1 3]);
z  = permute(z,[2 1 3]);
Rhoe = permute(Rhoe,[2 1 3]);
Rhoi = permute(Rhoi,[2 1 3]);
Bx = permute(Bx,[2 1 3]);
By = permute(By,[2 1 3]);
Bz = permute(Bz,[2 1 3]);
Ex = permute(Ex,[2 1 3]);
Ey = permute(Ey,[2 1 3]);
Ez = permute(Ez,[2 1 3]);
Uex = permute(Uex,[2 1 3]);
Uey = permute(Uey,[2 1 3]);
Uez = permute(Uez,[2 1 3]);
Uix = permute(Uix,[2 1 3]);
Uiy = permute(Uiy,[2 1 3]);
Uiz = permute(Uiz,[2 1 3]);

Coef = 1.0;
uxv= interp3(x, y, z, Uix, Coef*xq, Coef*yq, Coef*zq);
uyv= interp3(x, y, z, Uiy, Coef*xq, Coef*yq, Coef*zq);
uzv= interp3(x, y, z, Uiz, Coef*xq, Coef*yq, Coef*zq);

% Transform vectors into local coordinate system
uM = Inf(size(xq)); uL = uM; uN = uM;

% This could potentially be improved!
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      uM(ix,iy) = [uxv(ix,iy) uyv(ix,iy) uzv(ix,iy)]*dM(:,ix,iy);
      uL(ix,iy) = [uxv(ix,iy) uyv(ix,iy) uzv(ix,iy)]*dL(:,ix,iy);
      uN(ix,iy) = [uxv(ix,iy) uyv(ix,iy) uzv(ix,iy)]*dN(:,ix,iy);
   end
end


%
figure
a = surf(xq,yq,zq,uL);  axis equal
colorbar; colormap(jet);
a.LineStyle = 'none';
view([-67.9 14.8]);


%% movie
filename = '~/Ganymede/newPIC/run_G8_newPIC/3d_fluid_35.outs';
[~,~,fileinfo] = read_data(filename,'verbose',false);
npict = fileinfo.npictinfiles;

Coef = 1.00;
%height = 1;

% Create video
% v = VideoWriter('~/Ganymede/PC_ViSlice.avi');
% v.FrameRate = 10;
% v.open

% create new figure with specified size
hfig = figure(4);
set(hfig,'position', [10, 10, 800, 200]) 
colormap(jet);

for ipict = 1:1%npict
   [filehead,data] = read_data(filename,'verbose',false,'npict',ipict);
   
   x = data.file1.x(:,:,:,1);
   y = data.file1.x(:,:,:,2);
   z = data.file1.x(:,:,:,3);
   
%    uxe = data.file1.w(:,:,:,9);
%    uye = data.file1.w(:,:,:,10);
%    uze = data.file1.w(:,:,:,11);
   uxi = data.file1.w(:,:,:,12);
   uyi = data.file1.w(:,:,:,13);
   uzi = data.file1.w(:,:,:,14);
   pe  = data.file1.w(:,:,:,15);
   
   % From ndgrid to meshgrid format
   x  = permute(x,[2 1 3]);
   y  = permute(y,[2 1 3]);
   z  = permute(z,[2 1 3]);
%    uex = permute(uxe,[2 1 3]);
%    uey = permute(uye,[2 1 3]);
%    uez = permute(uze,[2 1 3]);
   uix = permute(uxi,[2 1 3]);
   uiy = permute(uyi,[2 1 3]);
   uiz = permute(uzi,[2 1 3]);   
   pe = permute(pe,[2 1 3]);
  
%    uxv= interp3(x, y, z, uxe, Coef*xq, Coef*yq, Coef*zq);
%    uyv= interp3(x, y, z, uye, Coef*xq, Coef*yq, Coef*zq); 
%    uzv= interp3(x, y, z, uze, Coef*xq, Coef*yq, Coef*zq); 
%    pev= interp3(x, y, z, pe, Coef*xq, Coef*yq, Coef*zq);
   
   uxv= interp3(x, y, z, uix, Coef*xq, Coef*yq, Coef*zq);
   uyv= interp3(x, y, z, uiy, Coef*xq, Coef*yq, Coef*zq); 
   uzv= interp3(x, y, z, uiz, Coef*xq, Coef*yq, Coef*zq); 
      
   % Transform vectors into local coordinate system
   u1 = Inf(size(xq)); u2 = u1; u3 = u1;
      
   % This could potentially be improved!
   for ix=1:size(xq,1)
      for iy=1:size(xq,2)
         uM(ix,iy) = [uxv(ix,iy) uyv(ix,iy) uzv(ix,iy)]*dM(:,ix,iy);
         uL(ix,iy) = [uxv(ix,iy) uyv(ix,iy) uzv(ix,iy)]*dL(:,ix,iy);
         uN(ix,iy) = [uxv(ix,iy) uyv(ix,iy) uzv(ix,iy)]*dN(:,ix,iy);
      end
   end
 
   subplot_tight(1,3,1);
   surf(xq,yq,zq,uM,'Linestyle','none'); colorbar;
   axis tight equal
   axis([min(x(:)) max(x(:)) ymin ymax zmin+0.1 zmax]);
   %caxis([150 1000])
   caxis([-150 150])
   xlabel('x [R_G]'); ylabel('y [R_G]'); zlabel('z [R_G]')
   title('U_M [km/s]')
   view([-67.9 14.8]);
   %h.Position(4) = height;
   
   subplot_tight(1,3,2);
   surf(xq,yq,zq,uL,'Linestyle','none'); colorbar;
   axis tight equal
   axis([min(x(:)) max(x(:)) ymin ymax zmin+0.1 zmax]);
   %caxis([-800 800])
   caxis([-200 200])
   xlabel('x [R_G]'); ylabel('y [R_G]'); zlabel('z [R_G]')
   title('U_L [km/s]')
   view([-67.9 14.8]);
      
   subplot_tight(1,3,3); 
   surf(xq,yq,zq,uN,'Linestyle','none'); colorbar;
   axis tight equal
   axis([min(x(:)) max(x(:)) ymin ymax zmin+0.1 zmax]);
   %caxis([0.2 5]);
   caxis([-80 100])
   xlabel('y [R_G]'); ylabel('z [R_G]'); zlabel('z [R_G]')
   %title('Pe [nPa]')
   title('U_N [km/s]')
   view([-67.9 14.8]);
   
   
   dim = [0.1 0.113 0.05 0.02];
   str = sprintf('it=%d, time=%.1fs',filehead.it,filehead.time);
   a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
      'FontWeight','bold','EdgeColor','none');
   
%    frame = getframe(gcf);
%    
%    set(gca,'nextplot','replacechildren');
%    
%    clf;
%    writeVideo(v,frame);
   
end

% v.close
% close(4)