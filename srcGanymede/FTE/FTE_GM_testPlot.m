

filename = '~/Ganymede/newPIC/test/box_Hall.out'; %GM symmetric test
filename = '~/Ganymede/newPIC/test/box_ideal.out'; %GM symmetric test

s = 1;
xThres = -1.5;
[x3bc,y3bc,z3bc] = find_bz0_boundary( filename,s,xThres );

ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

%% Generate mesh points from fitted surface and calculate LMN directions
ymin = -1.2; ymax = 1.2; zmin = -0.8+1/16; zmax = 0.8;
dy = 1/32; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

xq = fitresult(yq,zq);

% Calculate the normal direction to the fitted surface
[V, W] = differentiate(fitresult, yq, zq);
U = -ones(size(V));

% get the three local directions
% dipole-direction unit vector
unitDipole = [0 0 1];
% Initialize local vectors: d1-> M d2->L d3-> N
dL = Inf(3,size(xq,1),size(xq,2));
dM = dL; dN = dL;

% This part could potentially be optimized!
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      dN(:,ix,iy) = [U(ix,iy) V(ix,iy) W(ix,iy)];
      dM(:,ix,iy) = cross(dN(:,ix,iy),unitDipole);
      dL(:,ix,iy) = cross(dM(:,ix,iy),dN(:,ix,iy));
      
      % Normalization
      dL(:,ix,iy) = dL(:,ix,iy) / norm(dL(:,ix,iy));
      dM(:,ix,iy) = dM(:,ix,iy) / norm(dM(:,ix,iy));
      dN(:,ix,iy) = dN(:,ix,iy) / norm(dN(:,ix,iy));
   end
end


%%
% Offset by a distance of 1 cell (depends on grid resolution)
dn = dy;

Coef = 1.00;

% create new figure with specified size
hfig = figure;
set(hfig,'position', [10, 10, 800, 520]) 
colormap(jet);

[filehead,data] = read_data(filename,'verbose',false,'npict',1);

x = data.file1.x(:,:,:,1);
y = data.file1.x(:,:,:,2);
z = data.file1.x(:,:,:,3);

ux = data.file1.w(:,:,:,2);
uy = data.file1.w(:,:,:,3);
uz = data.file1.w(:,:,:,4);
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

% % Transform vectors into local coordinate system
% uL = Inf(size(xq)); uM = uL; uN = uL;
% bL = Inf(size(xq)); bM = bL; bN = bL;
% jL = Inf(size(xq)); jM = jL; jN = jL;
% 
% % This could potentially be improved!
% for ix=1:size(xq,1)
%    for iy=1:size(xq,2)
%       bL(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
%       bM(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
%       bN(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dN(:,ix,iy);
%       jL(ix,iy) = [jxv(ix,iy) jyv(ix,iy) jzv(ix,iy)]*dL(:,ix,iy);
%       jM(ix,iy) = [jxv(ix,iy) jyv(ix,iy) jzv(ix,iy)]*dM(:,ix,iy);
%       jN(ix,iy) = [jxv(ix,iy) jyv(ix,iy) jzv(ix,iy)]*dN(:,ix,iy);
%    end
% end

% Bx (positive pointing upstream)
subplot_tight(3,3,1);
contourf(yq,zq,bxv,50,'Linestyle','none'); colorbar;
axis tight equal
set(gca,'Xdir','reverse')
%caxis([-60 60])
xlabel('y [R_G]'); ylabel('z [R_G]');
title('B_x [nT]')

% By (positive pointing in roughly y direction)
subplot_tight(3,3,2);
contourf(yq,zq,byv,50,'Linestyle','none'); colorbar;
axis tight equal
set(gca,'Xdir','reverse')
%caxis([-60 30])
xlabel('y [R_G]'); ylabel('z [R_G]');
title('B_y [nT]')

% Bz
subplot_tight(3,3,3);
contourf(yq,zq,bzv,50,'Linestyle','none'); colorbar;
axis tight equal
set(gca,'Xdir','reverse')
%caxis([0.1 0.7])

xlabel('y [R_G]'); ylabel('z [R_G]');
title('B_z [nT]')


% jx
subplot_tight(3,3,4);
contourf(yq,zq,jxv,50,'Linestyle','none'); colorbar;
axis tight equal
set(gca,'Xdir','reverse')
%caxis([-200 200])
xlabel('y [R_G]'); ylabel('z [R_G]');
title('j_x [\mu A/m^2]')

% jy
subplot_tight(3,3,5);
contourf(yq,zq,jyv,50,'Linestyle','none'); colorbar;
axis tight equal
set(gca,'Xdir','reverse')
%caxis([-100 100])
xlabel('y [R_G]'); ylabel('z [R_G]');
title('j_y [\mu A/m^2]')


% jz
subplot_tight(3,3,6);
contourf(yq,zq,jzv,'Linestyle','none'); colorbar;
axis tight equal
set(gca,'Xdir','reverse')
%caxis([0.2 4]);
%caxis([2 12]);
xlabel('y [R_G]'); ylabel('z [R_G]');
title('j_z [\mu A/m^2]')



% ux
subplot_tight(3,3,7);
contourf(yq,zq,uxv,50,'Linestyle','none'); colorbar;
axis tight equal
set(gca,'Xdir','reverse')
%caxis([-200 200])
xlabel('y [R_G]'); ylabel('z [R_G]');
title('u_x [km/s]')

% uy
subplot_tight(3,3,8);
contourf(yq,zq,uyv,50,'Linestyle','none'); colorbar;
axis tight equal
set(gca,'Xdir','reverse')
%caxis([-100 100])
xlabel('y [R_G]'); ylabel('z [R_G]');
title('u_y [km/s]')


% uz
subplot_tight(3,3,9);
contourf(yq,zq,uzv,'Linestyle','none'); colorbar;
axis tight equal
set(gca,'Xdir','reverse')
%caxis([0.2 4]);
%caxis([2 12]);
xlabel('y [R_G]'); ylabel('z [R_G]');
title('u_z [km/s]')


dim = [0.17 0.05 0.05 0.02];
str = sprintf('it=%d, time=%.1fs',...
   filehead.it,filehead.time);
a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
   'FontWeight','bold','EdgeColor','none');