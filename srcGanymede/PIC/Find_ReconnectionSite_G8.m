% Find reconnection site using div V in tangential directions
%
% Hongyang Zhou, hyzhou@umich.edu 03/12/2018

clear; clc
%%
filenamePC='~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600s.outs';
%filenamePC='~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600to1200.outs';
[~,~,fileinfo] = read_data(filenamePC,'verbose',false);
npict = fileinfo.npictinfiles;


% Take divergence only for the tangential components
%%   
% First, I need to fit the magnetopause surface
% filename = '~/Ganymede/newPIC/run_G8_newPIC/3d_G8_steady.outs'; % 3d GM outputs
% s = 0.8; % compact boundary factor [0,1]
% 
% [x3bc,y3bc,z3bc] = find_boundary_points( filename,s );

filename = '~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600s.outs'; %PC
%filename='~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600to1200.outs';
s = 1; % compact boundary factor [0,1]
xThres = -1.5;
[x3bc,y3bc,z3bc] = find_bz0_boundary( filename,s,xThres );

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

%ymin = -1.2; ymax = 1.2; zmin = -0.45; zmax = 0.7;
ymin = -1.2; ymax = 1.2; zmin = -0.8+1/16; zmax = 0.8;
dy = 1/32; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

xq = fitresult(yq,zq);

% Boundary normal coordinate
% dipole-direction unit vector
unitDipole = [18 -51.82 716.8]/sqrt(18^2+51.82^2+716.8^2);
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

%%
% v = VideoWriter('~/Ganymede/divVuL_G8_600to1200.avi');
% v.FrameRate = 10;
% v.open

% Offset by a distance of 1 cell (depends on grid resolution)
dn = dy;
e = 1.6022e-19; %[C]

fig = figure(4);
fig.Position = [200,200,800,500];
set(fig,'nextplot','replacechildren');
colormap(jet);
box on

ReconnectionCount = zeros(size(xq));

for ipict = 1:1%npict
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
   bxv= interp3(x, y, z, Bx, xq, yq, zq);
   byv= interp3(x, y, z, By, xq, yq, zq);
   bzv= interp3(x, y, z, Bz, xq, yq, zq);
   
   bL = sum([bxv(:) byv(:) bzv(:)]'.*reshape(dL,[3,numel(xq)]));
   bL = reshape(bL,size(xq));
   
   offset = zeros(size(xq));
   % Loop over the outer surfaces
   for isurf = 1:10
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
   
   %ueL = medfilt2(ueL);
   %ueL = ordfilt2(ueL,16,ones(4,4));
   
   % Doing gradient using chain rule
   [pz,py] = gradient(ueL,dy);
   duLdL = py .* squeeze(dL(2,:,:)) + pz .* squeeze(dL(3,:,:));
   
   [pz,py] = gradient(ueM,dy);
   duMdM = py .* squeeze(dM(2,:,:)) + pz .* squeeze(dM(3,:,:));
   
   divV = duLdL + duMdM;
   
   % Pick a threshold
   divVthres = 2.5e3;
   %k = find(abs(duLdL) > divVthres & abs(ueL)<50);
   k = find(divV > divVthres & abs(ueL)<50);
   
   ReconnectionCount(k) = ReconnectionCount(k) + 1;
   
   surf(xqNew,yqNew,zqNew,ueL,'Linestyle','none'); %c = colorbar;
   %c.Limits = [-650 650]; c.Ticks = [-600 -400 -200 0 200 400 600];
   colorbar; caxis([-650 650]);
   axis tight equal
   view([-78.3 11.6]); hold on
   scatter3(xqNew(k),yqNew(k),zqNew(k),'k+');
   
   xlabel('x [R_G]'); ylabel('y [R_G]'); zlabel('z [R_G]');
   title('UeL [km/s] and reconnection sites')
   
   xlim([-2.0 -1.1]);
   ylim([ymin ymax]); zlim([zmin zmax]);
   
   dim = [0.2 0.113 0.05 0.02];
   str = sprintf('it=%d, time=%.1fs',...
      filehead.it,filehead.time);
   a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
      'FontWeight','bold','EdgeColor','none');
%    
%    frame = getframe(gcf);
%    
%    clf;
%    writeVideo(v,frame);

end

% figure;
% contourf(yq,zq,ReconnectionCount)
% save('ReconnectionCount_G8_600to1200.mat','ReconnectionCount')
% 
% v.close
% close(4)

% Find the outermost points in each y slices
OuterIndex = zeros(size(ueL,1),1);
%OuterPt = OuterIndex;

for ix = 1:size(ueL,1)
   [~,OuterIndex(ix)] = min(xqNew(ix,:));
   %OuterPt(ix) = xqNew(ix,OuterIndex(ix));
   xqOut(ix) = xqNew(ix,OuterIndex(ix));
   yqOut(ix) = yqNew(ix,OuterIndex(ix));
   zqOut(ix) = zqNew(ix,OuterIndex(ix));
end

scatter3(xqOut(:),yqOut(:),zqOut(:),'+')