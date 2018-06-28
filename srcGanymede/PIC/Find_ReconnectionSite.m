% Find reconnection site inside PIC box.
%
% Methods:
% 1. Find the large Eprime spots;
% 2. Find the large \vec{E}^prime \cdot \vec{j} spots;
% 3. Find the large J/B spots;
% 4. Find the large \div \cdot\vec{v} spots?
%
% Due to the existence of LHDI, do smoothing before.
%
% Note that some methods require steady state 3d GM outputs.
%
% Hongyang Zhou, hyzhou@umich.edu 02/22/2018

clear; clc
%% Parameters
e = 1.6022e-19; %[C]

method = 9;

ipict = 256;

%% Read data from PC
%filenamePC='~/Ganymede/newPIC/run_G8_newPIC/3d_fluid_35.outs';
filenamePC='~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600s.outs';
[filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);

data = data.file1;
x = data.x(:,:,:,1);       % [Rg]
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);
ne = data.w(:,:,:,1)*1e6;    % [/m^3]
ni = data.w(:,:,:,2)*1e6/14; % [/m^3] 
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
Pe = data.w(:,:,:,15);     % [nPa]
Pi = data.w(:,:,:,16);

% The original data is saved in ndgrid format. For streamline and
% isonormals functions, the input should be in meshgrid format.
x  = permute(x,[2 1 3]);
y  = permute(y,[2 1 3]);
z  = permute(z,[2 1 3]);
ne = permute(ne,[2 1 3]);
ni = permute(ni,[2 1 3]);
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
Pe  = permute(Pe,[2 1 3]);
Pi  = permute(Pi,[2 1 3]);

% Electron pressure gradient + electron inertia term [mV/m]
UeCrossB = cross([Uex(:) Uey(:) Uez(:)],[Bx(:) By(:) Bz(:)])*1e-3;

% Current density
J = e*(ni(:).*[Uix(:) Uiy(:) Uiz(:)] - ...
   ne(:).*[Uex(:) Uey(:) Uez(:)])*1e9; %[\mu A/m^2]

%% 
switch method
   case 1    
      Eprime = [Ex(:) Ey(:) Ez(:)] + UeCrossB;
      Eprime = sqrt(sum(Eprime.^2,2));
      Eprime = reshape(Eprime,size(Ex));
      
      %
      % This distribution plot is really interesting
      % figure;
      % histogram(Eprime(:));
      
      % Pick a threshold
      EprimeThres = 40; %[mV/m]
      k = find(Eprime > EprimeThres & Eprime < 41);

   case 2
      Eprime = [Ex(:) Ey(:) Ez(:)] + UeCrossB;
      
      % Dissipation term
      Dissipation = dot(Eprime',J');
      
      %figure;
      %histogram(Dissipation(:));
      
      % Pick a threshold
      % mean(Dissipation(:))+9*var(Dissipation(:))
      DissipationThres = 3.0; %[]
      k = find(Dissipation > DissipationThres);
      
   case 3
      J = sqrt(sum(J.^2,2));
      B = sqrt(Bx.^2 + By.^2 + Bz.^2);
      JoverB = J./B(:);
      %figure;
      %histogram(JoverB);
      
      % Pick a threshold
      JBthres = 4e-3;
      k = find(JoverB > JBthres);
   
   case 4
      divV = divergence(x,y,z,Uex,Uey,Uez);
      
      figure
      histogram(divV)
      
      % Pick a threshold
      divVthres = -6000;
      k = find(divV < divVthres);
   case 5
      divV = divergence(x,y,z,Uix,Uiy,Uiz);
      
      figure
      histogram(divV)
      
      % Pick a threshold
      divVthres = -1000;
      k = find(divV < divVthres);
      
   case 6
      % Take divergence only for the tangential components
      
      % First, I need to fit the magnetopause surface
      filename = '~/Ganymede/newPIC/run_G8_newPIC/3d_G8_steady.outs'; % 3d GM outputs
      s = 0.5; % compact boundary factor [0,1]

      [x3bc,y3bc,z3bc] = find_boundary_points( filename,s );
      
      % Set up fittype and options.
      ft = fittype( 'poly55' );
      
      % Fit model to data.
      [fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );
      
      ymin = -1.5; ymax = 1.5; zmin = -1; zmax = 1;
      dy = 1/30; dz = dy;
      %[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);
      [yq,zq] = meshgrid(ymin:dy:ymax,zmin:dz:zmax);
      
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
  
      % Interpolate onto the magnetopause surface
      Coef = 1.0;
      uxv= interp3(x, y, z, Uix, Coef*xq, Coef*yq, Coef*zq);
      uyv= interp3(x, y, z, Uiy, Coef*xq, Coef*yq, Coef*zq);
      uzv= interp3(x, y, z, Uiz, Coef*xq, Coef*yq, Coef*zq);
      
      % Transform vectors into local coordinate system
      uM = Inf(size(xq)); uL = uM; uN = uM;
      
      %
      for ix=1:size(xq,1)
         for iy=1:size(xq,2)
            uM(ix,iy) = [uxv(ix,iy) uyv(ix,iy) uzv(ix,iy)]*dM(:,ix,iy);
            uL(ix,iy) = [uxv(ix,iy) uyv(ix,iy) uzv(ix,iy)]*dL(:,ix,iy);
            uN(ix,iy) = [uxv(ix,iy) uyv(ix,iy) uzv(ix,iy)]*dN(:,ix,iy);
         end
      end
      
      % Doing gradient using chain rule
      [pz,py] = gradient(uL,dy);
      duLdL = py .* squeeze(dL(2,:,:)) + pz .* squeeze(dL(3,:,:));
      
      [pz,py] = gradient(uM,dy);
      duMdM = py .* squeeze(dM(2,:,:)) + pz .* squeeze(dM(3,:,:));
      
      divV = duLdL + duMdM;
      
      % Pick a threshold
      divVthres = 300;
      k = find(duLdL > divVthres & xq<-1.6);
      %k = find(divV > divVthres & xq<-1.6);      
      
%       figure;
%       surf(xq,yq,zq,divV,'Linestyle','none'); colorbar;
%       %axis tight equal
%       
%       figure;
%       surf(xq,yq,zq,duLdL,'Linestyle','none'); colorbar;
% 
%       figure;
%       surf(xq,yq,zq,duMdM,'Linestyle','none'); colorbar;   
%       
      figure;
      surf(xq,yq,zq,uL,'Linestyle','none'); colorbar; hold on
      axis tight equal
      scatter3(xq(k),yq(k),zq(k),'r');
      
      %me       = 9.1094e-31;   %[kg]
      mi = 14; % [amu]
      mp       = 1.6726*1e-27; %[kg]
      me       = mp * mi/100;
      e        = 1.6022*1e-19; %[C]
      B = sqrt(Bx.^2 + By.^2 + Bz.^2);
      Ue = sqrt(Uex.^2 + Uey.^2 + Uez.^2);
      Ui = sqrt(Uix.^2 + Uiy.^2 + Uiz.^2);
      Re = me.*Ue*1e3./(e*B*1e-9);
      Re= interp3(x, y, z, Re, Coef*xq, Coef*yq, Coef*zq);
      Ri = mi*mp.*Ui*1e3./(e*B*1e-9);
      Ri= interp3(x, y, z, Ri, Coef*xq, Coef*yq, Coef*zq);      
      figure;
      surf(xq,yq,zq,2*pi/0.85*sqrt(Re.*Ri)/2634000,'Linestyle','none'); colorbar;
      
      return
   case 7
      % \partial B_N /\partial L, \partial B_L /\partial N
      % First, I need to fit the magnetopause surface
      filename = '~/Ganymede/newPIC/run_G8_newPIC/3d_G8_steady.outs'; % 3d GM outputs
      s = 0.5; % compact boundary factor [0,1]
      
      [x3bc,y3bc,z3bc] = find_boundary_points( filename,s );
      
      % Set up fittype and options.
      ft = fittype( 'poly55' );
      
      % Fit model to data.
      [fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );
      
      ymin = -1.5; ymax = 1.5; zmin = -1; zmax = 1;
      dy = 1/30; dz = dy;
      %[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);
      [yq,zq] = meshgrid(ymin:dy:ymax,zmin:dz:zmax);
      
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
      
      % Interpolate onto the magnetopause surface
      Coef = 1.0;
      bxv= interp3(x, y, z, Bx, Coef*xq, Coef*yq, Coef*zq);
      byv= interp3(x, y, z, By, Coef*xq, Coef*yq, Coef*zq);
      bzv= interp3(x, y, z, Bz, Coef*xq, Coef*yq, Coef*zq);
      
      % Transform vectors into local coordinate system
      bM = Inf(size(xq)); bL = bM; bN = bM;
      
      %
      for ix=1:size(xq,1)
         for iy=1:size(xq,2)
            bM(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
            bL(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
            bN(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dN(:,ix,iy);
         end
      end
      
      % Doing gradient using chain rule
      [pz,py] = gradient(bN,dy);
      dbNdL = py .* squeeze(dL(2,:,:)) + pz .* squeeze(dL(3,:,:));
      
      figure;
      surf(xq,yq,zq,dbNdL,'Linestyle','none'); colorbar;
      
      return
   case 8
      % E_perp
      Eperp = cross(J,[Bx(:) By(:) Bz(:)])./ (ni(:)*e)*1e-12; %[mV/m]
      Eperp = sqrt(sum(Eperp.^2,2));      
      Eperp = reshape(Eperp,size(Ex));
      
      sliceomatic(Eperp)
      return
   case 9
      % Ion [mV/m]
      UiCrossB = cross([Uix(:) Uiy(:) Uiz(:)],[Bx(:) By(:) Bz(:)])*1e-3;
      
      Eprime = [Ex(:) Ey(:) Ez(:)] + UiCrossB;
      Eprime = sqrt(sum(Eprime.^2,2));
      Eprime = reshape(Eprime,size(Ex));
      
      %
      % This distribution plot is really interesting
      % figure;
      % histogram(Eprime(:));
      
      % Pick a threshold
      EprimeThres = 90; %[mV/m]
      k = find(Eprime > EprimeThres);      
      
end

figure;
scatter3(x(k),y(k),z(k),'*')
axis equal
axis([min(x(:)) max(x(:)) min(y(:)) max(y(:)) min(z(:)) max(z(:))]);
view([-78.3 11.6]);
