% Reconnection site bin plot for G28
%
% Hongyang Zhou, hyzhou@umich.edu 03/06/2018

%%
clear;clc
% Find reconnection site using div V in tangential directions

filenamePC='~/Ganymede/newPIC/G28_PIC_theta51/3d_fluid_600s.outs';
[~,~,fileinfo] = read_data(filenamePC,'verbose',false);
npict = fileinfo.npictinfiles;

% Take divergence only for the tangential components
   
% First, I need to fit the magnetopause surface
% filename = '~/Ganymede/newPIC/run_G8_newPIC/3d_G8_steady.outs'; % 3d GM outputs
% s = 1.0; % compact boundary factor [0,1]
% 
% [x3bc,y3bc,z3bc] = find_boundary_points( filename,s );

filename = '~/Ganymede/newPIC/G28_PIC_theta51/3d_fluid_600s.outs'; %PC
s = 1; % compact boundary factor [0,1]
xThres = -1.5;
[x3bc,y3bc,z3bc] = find_bz0_boundary( filename,s,xThres );

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

% ymin = -1.5; ymax = 1.5; zmin = -1; zmax = 1;
ymin = -1.2; ymax = 1.2; zmin = -0.8+1/16; zmax = 0.8;
dy = 1/32; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

xq = fitresult(yq,zq);

% Boundary normal coordinate
% dipole-direction unit vector
unitDipole = [19.26 -16.54 716.8]/sqrt(19.26^2+16.54^2+716.8^2);
% Initialize local vectors
dM = Inf(3,size(xq,1),size(xq,2)); dL = dM; dN = dM;

% Calculate the normal direction to the fitted surface
[V, W] = differentiate(fitresult, yq, zq);
U = -ones(size(V));

% Get the three local directions
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      dN(:,ix,iy) = [U(ix,iy) V(ix,iy) W(ix,iy)];
      dM(:,ix,iy) = cross(dN(:,ix,iy),unitDipole);
      dL(:,ix,iy) = cross(dM(:,ix,iy),dN(:,ix,iy));
      
      % Normalization
      dM(:,ix,iy) = dM(:,ix,iy) / norm(dM(:,ix,iy));
      dL(:,ix,iy) = dL(:,ix,iy) / norm(dL(:,ix,iy));
      dN(:,ix,iy) = dN(:,ix,iy) / norm(dN(:,ix,iy));
   end
end

ReconnectionCount = zeros(size(xq));

% Offset by a distance of 1 cell (depends on grid resolution)
dn = dy;

for ipict = 1:npict
   fprintf('ipict=%d\n',ipict);
   [filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);
   
   data = data.file1;
   x = data.x(:,:,:,1);       % [Rg]
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);
   Bx = data.w(:,:,:,3);      % [nT]
   By = data.w(:,:,:,4);
   Bz = data.w(:,:,:,5);
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
   Bx = permute(Bx,[2 1 3]);
   By = permute(By,[2 1 3]);
   Bz = permute(Bz,[2 1 3]);
   Uex = permute(Uex,[2 1 3]);
   Uey = permute(Uey,[2 1 3]);
   Uez = permute(Uez,[2 1 3]);
   Uix = permute(Uix,[2 1 3]);
   Uiy = permute(Uiy,[2 1 3]);
   Uiz = permute(Uiz,[2 1 3]);
   
   
   % Get the starting surface points value
   isurf = -5;
   bxv= interp3(x, y, z, Bx, xq + isurf*dn*squeeze(dN(1,:,:)),...
      yq + isurf*dn*squeeze(dN(2,:,:)),...
      zq + isurf*dn*squeeze(dN(3,:,:)));
   byv= interp3(x, y, z, By, xq + isurf*dn*squeeze(dN(1,:,:)),...
      yq + isurf*dn*squeeze(dN(2,:,:)), ...
      zq + isurf*dn*squeeze(dN(3,:,:)));
   bzv= interp3(x, y, z, Bz, xq + isurf*dn*squeeze(dN(1,:,:)),...
      yq + isurf*dn*squeeze(dN(2,:,:)), ...
      zq + isurf*dn*squeeze(dN(3,:,:)));
   
   bL = sum([bxv(:) byv(:) bzv(:)]'.*reshape(dL,[3,numel(xq)]));
   bL = reshape(bL,size(xq));
   
   offset = zeros(size(xq));
   % Loop over the outer surfaces
   for isurf = -4:10
      bxv= interp3(x, y, z, Bx, xq + isurf*dn*squeeze(dN(1,:,:)),...
         yq + isurf*dn*squeeze(dN(2,:,:)),...
         zq + isurf*dn*squeeze(dN(3,:,:)));
      byv= interp3(x, y, z, By, xq + isurf*dn*squeeze(dN(1,:,:)),...
         yq + isurf*dn*squeeze(dN(2,:,:)), ...
         zq + isurf*dn*squeeze(dN(3,:,:)));
      bzv= interp3(x, y, z, Bz, xq + isurf*dn*squeeze(dN(1,:,:)),...
         yq + isurf*dn*squeeze(dN(2,:,:)), ...
         zq + isurf*dn*squeeze(dN(3,:,:)));
      
      % Get bL on the outside surface
      bLOut = sum([bxv(:) byv(:) bzv(:)]'.*reshape(dL,[3,numel(xq)]));
      bLOut = reshape(bLOut,size(xq));

      % Assuming bL changes monotonically     
      offset = offset + dn * (bL.*bLOut<0) .* (isurf-1+bL./(bL-bLOut));     
      bL = bLOut;
   end
   
   offset(offset==0) = nan; % Remove the unecessary points
   
   xqNew = xq + offset.*squeeze(dN(1,:,:));
   yqNew = yq + offset.*squeeze(dN(2,:,:));
   zqNew = zq + offset.*squeeze(dN(3,:,:));
   
   Uexv= interp3(x, y, z, Uex, xqNew, yqNew, zqNew);
   Ueyv= interp3(x, y, z, Uey, xqNew, yqNew, zqNew); 
   Uezv= interp3(x, y, z, Uez, xqNew, yqNew, zqNew);
   
   ueL = Inf(size(xq)); ueM = ueL; %ueN = ueL;
   
   % This could potentially be improved!
   for ix=1:size(xq,1)
      for iy=1:size(xq,2)
         ueL(ix,iy) = [Uexv(ix,iy) Ueyv(ix,iy) Uezv(ix,iy)]*dL(:,ix,iy);
         ueM(ix,iy) = [Uexv(ix,iy) Ueyv(ix,iy) Uezv(ix,iy)]*dM(:,ix,iy);
         %ueN(ix,iy) = [Uexv(ix,iy) Ueyv(ix,iy) Uezv(ix,iy)]*dN(:,ix,iy);
      end
   end   
   
   % Doing gradient using chain rule
   [pz,py] = gradient(ueL,dy);
   duLdL = py .* squeeze(dL(2,:,:)) + pz .* squeeze(dL(3,:,:));
   
   [pz,py] = gradient(ueM,dy);
   duMdM = py .* squeeze(dM(2,:,:)) + pz .* squeeze(dM(3,:,:));
   
   divV = duLdL + duMdM;
   
   % Pick a threshold
   divVthres = 2.5e3;
   k = find(divV > divVthres & abs(ueL)<50);
   
   ReconnectionCount(k) = ReconnectionCount(k) + 1;
end

figure;
contourf(yq,zq,ReconnectionCount)
save('ReconnectionCount_G28.mat','ReconnectionCount')