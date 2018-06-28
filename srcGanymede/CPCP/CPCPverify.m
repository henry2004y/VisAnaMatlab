% Script for 3d box outputs.
%
% Hongyang Zhou, hyzhou@umich.edu 01/18/2018

clear; clc
%%
Rg = 2634000; %[m], radius of Ganymede
e  = 1.60217662e-19; % [C], electron charge

filename= '~/Ganymede/MOP2018/runG28_PIC_CFL04_50s_test/GM/box_var_6_n50000_79780.outs';

[filehead,data] = read_data(filename,'verbose',false,'npict',51);

data = data.file1;

x = data.x(:,:,:,1);
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);
status = data.w(:,:,:,14);

% Topology plots
figure;
contourf(x(:,:,1),y(:,:,1),status(:,:,1),'EdgeColor','none'); colorbar

figure;
contourf(x(:,:,2),y(:,:,2),status(:,:,2),'EdgeColor','none'); colorbar

figure;
contourf(x(:,:,3),y(:,:,3),status(:,:,3),'EdgeColor','none'); colorbar

%% Check the electric field integral for 3 cuts

for icut = 1:3
   
   % Get E field: E = -uxB
   rho= data.w(:,:,icut,1);rho = rho(:);
   ux = data.w(:,:,icut,2); ux = ux(:);
   uy = data.w(:,:,icut,3); uy = uy(:);
   uz = data.w(:,:,icut,4); uz = uz(:);
   
   u = [ux uy uz];
   
   Bx = data.w(:,:,icut,5); Bx = Bx(:);
   By = data.w(:,:,icut,6); By = By(:);
   Bz = data.w(:,:,icut,7); Bz = Bz(:);
   
   B = [Bx By Bz];
   
   jx = data.w(:,:,icut,11); jx = jx(:);
   jy = data.w(:,:,icut,12); jy = jy(:);
   jz = data.w(:,:,icut,13); jz = jz(:);
   
   J = [jx jy jz];
   
   % Calculate electron bulk velocity in Hall MHD
   % Transform the hall velocity into [km/s]
   uex = ux - jx./rho * 1e-6 /e * 1e-3 * 1e-6;
   uey = uy - jy./rho * 1e-6 /e * 1e-3 * 1e-6;
   uez = uz - jz./rho * 1e-6 /e * 1e-3 * 1e-6;
   
   ue = [uex uey uez];
   
   E = -cross(ue,B);
   E = reshape(E,size(data.w,1),size(data.w,2),3);
   
   Ex = E(:,:,1);
   Ey = E(:,:,2);
   Ez = E(:,:,3);
   
   statusIn = (status(:,:,icut)==2);
   
   xIn = x(:,:,icut);
   yIn = y(:,:,icut);
   
   xIn = xIn(statusIn);
   yIn = yIn(statusIn);
   
   % Find the boundary of magnetosphere
   k = boundary(xIn,yIn,1);
   %figure(1); hold on
   %plot(xIn(k),yIn(k));
   
   % Integrate along the boundary
   xIn = xIn(k);
   yIn = yIn(k);
   Ex = Ex(statusIn);
   Ex = Ex(k);
   Ey = Ey(statusIn);
   Ey = Ey(k);
   
   % Calculate cross polar cap potential (SI units!)
   % Let the starting point have reference absolute potential 0
   Ediff = Inf(numel(xIn)-1,1);
   for iE=1:numel(xIn)-1
      Ediff(iE) = ( (Ex(iE)+Ex(iE+1))/2*(xIn(iE+1)-xIn(iE)) + ...
         (Ey(iE)+Ey(iE+1))/2*(yIn(iE+1)-yIn(iE)) )*Rg*1e3*1e-9;
   end
   
   EPotential = cumsum([0;Ediff]);
   
   figure(4); hold on
   plot(EPotential)

end
