% Find LMN direction using Minimum Derivative Direction method.

clear; clc;% close all
%%

filename = '~/Ganymede/MOP2018/runG8_PIC_1200s/3d_t=280.out';
s = 0.9; % compact boundary factor [0,1]

% PC upstream box outputs
PCdir = '~/Ganymede/MOP2018/runG8_PIC_1200s/PC';
PCfile= '3d_var_region0_0_t00000135_n00002849.out';
%PCfile= '3d_var_region0_0_t00000444_n00008549.out';

% %% Find boundary points from steady state solution
% 
% [x3bc,y3bc,z3bc] = find_boundary_points( filename,s );
% 
% %% Fit the closed field line boundary with hypersurface
% 
% [fitresult,gof] = surface_fit(x3bc,y3bc,z3bc);
% 
% %% Generate mesh points from fitted surface and calculate LMN directions
% ymin = -1.2; ymax = 1.2; zmin = -0.6; zmax = 0.6;
% dy = 1/32; dz = dy;
% [yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);
% 
% xq = fitresult(yq,zq);
% 
% % Calculate the normal direction to the fitted surface
% [V, W] = differentiate(fitresult, yq, zq);
% U = -ones(size(V));
% 
% % get the three local directions
% % dipole-direction unit vector
% unitDipole = [18 -51.82 716.8]/sqrt(18^2+51.82^2+716.8^2);
% % Initialize local vectors: d1-> M d2->L d3-> N
% dL = Inf(3,size(xq,1),size(xq,2));
% dM = dL; dN = dL;
% 
% % This part could potentially be optimized!
% for ix=1:size(xq,1)
%    for iy=1:size(xq,2)
%       dN(:,ix,iy) = [U(ix,iy) V(ix,iy) W(ix,iy)];
%       dM(:,ix,iy) = cross(dN(:,ix,iy),unitDipole);
%       dL(:,ix,iy) = cross(dM(:,ix,iy),dN(:,ix,iy));
%       
%       % Normalization
%       dL(:,ix,iy) = dL(:,ix,iy) / norm(dL(:,ix,iy));
%       dM(:,ix,iy) = dM(:,ix,iy) / norm(dM(:,ix,iy));
%       dN(:,ix,iy) = dN(:,ix,iy) / norm(dN(:,ix,iy));
%    end
% end

%% Get 3D B field

[filehead,data] = read_data(fullfile(PCdir,PCfile),...
   'verbose',false,'npict',1);

bx_ = strcmpi('bx',filehead.wnames);
by_ = strcmpi('by',filehead.wnames);
bz_ = strcmpi('bz',filehead.wnames);

data = data.file1;
x = data.x(:,:,:,1);
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

bx = data.w(:,:,:,bx_);      % [nT]
by = data.w(:,:,:,by_);
bz = data.w(:,:,:,bz_);

% From ndgrid to meshgrid format
% x  = permute(x,[2 1 3]);
% y  = permute(y,[2 1 3]);
% z  = permute(z,[2 1 3]);
% bx = permute(bx,[2 1 3]);
% by = permute(by,[2 1 3]);
% bz = permute(bz,[2 1 3]);

% dx = x(1,3,1) - x(1,2,1);
% dy = y(3,1,1) - y(2,1,1);
% dz = z(1,1,3) - z(1,1,2);

dx = x(3,1,1) - x(2,1,1);
dy = y(1,3,1) - y(1,2,1);
dz = z(1,1,3) - z(1,1,2);


% Cut the matrix for testting
%x = x()

[dBxdx,dBxdy,dBxdz] = gradient(bx,dx,dy,dz);
[dBydx,dBydy,dBydz] = gradient(by,dx,dy,dz);
[dBzdx,dBzdy,dBzdz] = gradient(bz,dx,dy,dz);

matB = zeros(3,3,numel(dBxdx));

nI = size(dBxdx,1); nJ = size(dBxdx,2); nK = size(dBxdx,3);

for k=1:nK
   for j=1:nJ
      for i=1:nI
         iB = (k-1)*nI*nJ + (j-1)*nI + i;   
         matB(:,:,iB) = [dBxdx(i,j,k) dBydx(i,j,k) dBzdx(i,j,k);
            dBxdy(i,j,k) dBydy(i,j,k) dBzdy(i,j,k);
            dBxdz(i,j,k) dBydz(i,j,k) dBzdz(i,j,k)];
         matB(:,:,iB) = matB(:,:,iB)*matB(:,:,iB)';
      end
   end
end

%%
% Find the magnetopause locations
sig = sign(bz);
fx  = gradient(sig);
BCindex_ = find(fx~=0);

x3bc = x(BCindex_);
y3bc = y(BCindex_);
z3bc = z(BCindex_);
%matB = matB(:,:,BCindex_);

BCindex_ = z3bc<0.8 & z3bc>-0.75;
x3bc = x3bc(BCindex_);
y3bc = y3bc(BCindex_);
z3bc = z3bc(BCindex_);
%matB = matB(:,:,BCindex_);

figure(1)
scatter3(x3bc,y3bc,z3bc,'.'); axis equal

% x = F(y,z)
F = scatteredInterpolant(y3bc,z3bc,x3bc);

ymin = -1.2; ymax = 1.2; zmin = -0.6; zmax = 0.6;
dy = 1/30; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

xq = F(yq,zq);

figure(2)
surf(xq,yq,zq,'FaceAlpha',0.7);
axis equal; grid on; box on
xlabel('x [R_G]'); ylabel('y [R_G]'); zlabel('z [R_G]');
set(gca,'LineWidth',1.2,'FontSize',14)


%%

% For interpolation, I need ndgrid format!
matB_new = Inf([3,3,size(xq)]);

matB = reshape(matB,[3,3,size(x)]);
% 3D interpolant
for j=1:3; for i=1:3
   F = griddedInterpolant(x,y,z,squeeze(matB(i,j,:,:,:)),'cubic');
   matB_new(i,j,:,:,:) = F(xq,yq,zq);
end; end


%scatter3(x(:),y(:),z(:),'.')

% range_x = linspace(min(x(:)),max(x(:)),size(x,1));
% range_y = linspace(min(y(:)),max(y(:)),size(x,2));
% range_z = linspace(min(z(:)),max(z(:)),size(x,3));
%gv = {range_x,range_y,range_z};
%F = griddedInterpolant(gv,matB(1,1,:),'cubic');

%%
V = Inf(size(matB_new));
D = Inf(size(matB_new));

for j=1:size(matB_new,4); for i=1:size(matB_new,3)
   [V(:,:,i,j),D(:,:,i,j)] = eig(matB_new(:,:,i,j));
end; end

figure(3)
% quiver3(xq,yq,zq,...
%    squeeze(-V(1,3,:,:)),squeeze(-V(2,3,:,:)),squeeze(-V(3,3,:,:)),1)
quiver3(xq,yq,zq,...
   squeeze(V(1,3,:,:)),squeeze(V(2,3,:,:)),squeeze(V(3,3,:,:)),1)
axis equal
xlabel('x'); ylabel('y'); zlabel('z');

figure(4)
quiver3(xq,yq,zq,...
   squeeze(V(1,2,:,:)),squeeze(V(2,2,:,:)),squeeze(V(3,2,:,:)),1)
axis equal
xlabel('x'); ylabel('y'); zlabel('z');

figure(5)
quiver3(xq,yq,zq,...
   squeeze(V(1,1,:,:)),squeeze(V(2,1,:,:)),squeeze(V(3,1,:,:)),1)
axis equal
xlabel('x'); ylabel('y'); zlabel('z');

%%
return

V = Inf(3,3,size(matB,3)); D = Inf(3,3,size(matB,3));

for i=1:size(matB,3)
   [V(:,:,i),D(:,:,i)] = eig(matB(:,:,i));
end

figure(3)
for i=1:size(matB,3)
   quiver3(x3bc(i),y3bc(i),z3bc(i),-V(1,3,i),-V(2,3,i),-V(3,3,i),0.05,'k')
   hold on
end

xlabel('x'); ylabel('y'); zlabel('z');




