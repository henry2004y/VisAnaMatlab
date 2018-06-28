% Creating movies for combining clock angles and reconnection sites.
%
% Needs to be re-organized and fully tested.

clear; clc
%%
Coefin  = 0.92; % shrink factor
Coefout = 1.08; % expansion factor

%% GM

ipict = 1;
filenameGM = '~/Ganymede/newPIC/run_G8_newPIC/box_FTE_G8_1200s.outs';
[filehead,data] = read_data(filenameGM,'verbose',false,'npict',ipict);

data = data.file1;
x = data.x(:,:,:,1); 
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

VarIndex_ = strcmpi('status',filehead.wnames);
status = data.w(:,:,:,VarIndex_);

x3_1 = x(status==3 & x<0);
y3_1 = y(status==3 & x<0);
z3_1 = z(status==3 & x<0);

% Pick the compact boundary points coords.
s = 0.5; % compact boundary factor [0,1]
k31  = boundary(x3_1,y3_1,z3_1,s);
bc31 = unique(k31);

x3bc = x3_1(bc31);
y3bc = y3_1(bc31);
z3bc = z3_1(bc31);

rThres = 1.6; % [Rg] radius threshold for finding dayside boundary pts
xThres = -1.5;

% Find the outer boundary points
mapindex_ = (x3bc.^2 + y3bc.^2 + z3bc.^2) > rThres^2 & x3bc < xThres;
x3bc = x3bc(mapindex_);
y3bc = y3bc(mapindex_);
z3bc = z3bc(mapindex_);

%Fit the closed field line boundary with hypersurface

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

% Generate meshgrid from fitted surface
ymin = -1.; ymax = 1.; zmin = -0.5; zmax = 0.5;
dy = 1/30; dz = dy;
[yq,zq] = meshgrid(ymin:dy:ymax,zmin:dz:zmax);
xq = fitresult(yq,zq);

% Calculate the normal direction to the fitted surface
[V, W] = differentiate(fitresult, yq, zq);
U = -ones(size(V));

%% get the three local directions
% dipole-direction unit vector
unitDipole = [18 -51.82 716.8]/sqrt(18^2+51.82^2+716.8^2);
% Initialize local vectors
dM = Inf(3,size(xq,1),size(xq,2));
dL = dM; dN = dM;

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

%% Movie1
% Find reconnection site using Eprime

filenamePC='~/Ganymede/newPIC/run_G8_newPIC/3d_fluid_35.outs';
%filenamePC='~/Ganymede/newPIC/3d_theta51_test.outs';
[~,~,fileinfo] = read_data(filenamePC,'verbose',false);
npict = fileinfo.npictinfiles;

v = VideoWriter('~/Ganymede/Eprime.avi');
v.FrameRate = 10;
v.open

fig = figure(1); 
set(fig,'nextplot','replacechildren');
box on

for ipict = 1:npict
   [filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);
   
   data = data.file1;
   x = data.x(:,:,:,1);       % [Rg]
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);
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
   
   % Electron pressure gradient + electron intertia term [mV/m]
   UeCrossB = cross([Uex(:) Uey(:) Uez(:)],[Bx(:) By(:) Bz(:)])*1e-3;
   Eprime = [Ex(:) Ey(:) Ez(:)] + UeCrossB;
   
   Eprime = sqrt(sum(Eprime.^2,2));
   Eprime = reshape(Eprime,size(Ex));
   
   % Pick a threshold
   EprimeThres = 60; %[mV/m]
   k = find(Eprime > EprimeThres);
   
   scatter3(x(k),y(k),z(k),'.'); hold on
   axis equal
   axis([min(x(:)) max(x(:)) min(y(:)) max(y(:)) min(z(:)) max(z(:))]);
   view([-78.3 11.6]);
   
   dim = [0.2 0.113 0.05 0.02];
   str = sprintf('it=%d, time=%.1fs',...
      filehead.it,filehead.time);
   a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
      'FontWeight','bold','EdgeColor','none');
   
   frame = getframe(gcf);
   
   clf;
   writeVideo(v,frame);

end

v.close
close(1)

%% Movie2
% Find reconnection site using \vec{E}^prime\cdot\vec{J}

e = 1.6022e-19; %[C]

filenamePC='~/Ganymede/newPIC/run_G8_newPIC/3d_fluid_35.outs';
%filenamePC='~/Ganymede/newPIC/3d_theta51_test.outs';
[~,~,fileinfo] = read_data(filenamePC,'verbose',false);
npict = fileinfo.npictinfiles;

v = VideoWriter('~/Ganymede/EprimeJ.avi');
v.FrameRate = 10;
v.open

fig = figure(2); 
set(fig,'nextplot','replacechildren');
box on

for ipict = 1:npict
   [filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);
   
   data = data.file1;
   x = data.x(:,:,:,1);       % [Rg]
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);
   rhoe = data.w(:,:,:,1)*1e6;    % [/m^3]
   rhoi = data.w(:,:,:,2)*1e6/14; % [/m^3]
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
   rhoe = permute(rhoe,[2 1 3]);
   rhoi = permute(rhoi,[2 1 3]);
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
   
   % Electron pressure gradient + electron intertia term [mV/m]
   UeCrossB = cross([Uex(:) Uey(:) Uez(:)],[Bx(:) By(:) Bz(:)])*1e-3;
   
   % Current density
   J = e*(rhoi(:).*[Uix(:) Uiy(:) Uiz(:)] - ...
      rhoe(:).*[Uex(:) Uey(:) Uez(:)])*1e9; %[\mu A/m^2]
   
   Eprime = [Ex(:) Ey(:) Ez(:)] + UeCrossB;
   
   % Dissipation term
   Dissipation = dot(Eprime',J');
   
   
   % Pick a threshold
   %mean(Dissipation(:))+9*var(Dissipation(:))
   DissipationThres = 2.5; %[]
   k = find(Dissipation > DissipationThres);
       
   scatter3(x(k),y(k),z(k),'.'); hold on
   axis equal
   axis([min(x(:)) max(x(:)) min(y(:)) max(y(:)) min(z(:)) max(z(:))]);
   view([-78.3 11.6]);

   dim = [0.2 0.113 0.05 0.02];
   str = sprintf('it=%d, time=%.1fs',...
      filehead.it,filehead.time);
   a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
      'FontWeight','bold','EdgeColor','none');   
   
   frame = getframe(gcf);
   
   clf;
   writeVideo(v,frame);

end

v.close
close(2)


%% Movie3
% Find reconnection site using J/B

filenamePC='~/Ganymede/newPIC/run_G8_newPIC/3d_fluid_35.outs';
%filenamePC='~/Ganymede/newPIC/3d_theta51_test.outs';
[~,~,fileinfo] = read_data(filenamePC,'verbose',false);
npict = fileinfo.npictinfiles;

v = VideoWriter('~/Ganymede/JoverB.avi');
v.FrameRate = 10;
v.open

fig = figure(3); 
set(fig,'nextplot','replacechildren');
box on

for ipict = 1:npict
   [filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);
   
   data = data.file1;
   x = data.x(:,:,:,1);       % [Rg]
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);
   rhoe = data.w(:,:,:,1)*1e6;    % [/m^3]
   rhoi = data.w(:,:,:,2)*1e6/14; % [/m^3]
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
   rhoe = permute(rhoe,[2 1 3]);
   rhoi = permute(rhoi,[2 1 3]);
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
   
   % Electron pressure gradient + electron intertia term [mV/m]
   UeCrossB = cross([Uex(:) Uey(:) Uez(:)],[Bx(:) By(:) Bz(:)])*1e-3;
   
   % Current density
   J = e*(rhoi(:).*[Uix(:) Uiy(:) Uiz(:)] - ...
      rhoe(:).*[Uex(:) Uey(:) Uez(:)])*1e9; %[\mu A/m^2]
   
   J = sqrt(sum(J.^2,2));
   B = sqrt(Bx.^2 + By.^2 + Bz.^2);
   JoverB = J./B(:);
   
   % Pick a threshold
   JBthres = 4e-3;
   k = find(JoverB > JBthres);
      
   
   scatter3(x(k),y(k),z(k),'.'); hold on
   axis equal
   axis([min(x(:)) max(x(:)) min(y(:)) max(y(:)) min(z(:)) max(z(:))]);
   view([-78.3 11.6]);

   dim = [0.2 0.113 0.05 0.02];
   str = sprintf('it=%d, time=%.1fs',...
      filehead.it,filehead.time);
   a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
      'FontWeight','bold','EdgeColor','none');   
   
   frame = getframe(gcf);
   
   clf;
   writeVideo(v,frame);

end

v.close
close(3)


%% Movie4
% Find reconnection site using div V in tangential directions

filenamePC='~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600s.outs';
[~,~,fileinfo] = read_data(filenamePC,'verbose',false);
npict = fileinfo.npictinfiles;

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

v = VideoWriter('~/Ganymede/divVuL_new.avi');
v.FrameRate = 10;
v.open

fig = figure(4);
fig.Position = [200,200,800,500];
set(fig,'nextplot','replacechildren');
box on

for ipict = 538:538%npict
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
   %divV = pz .* squeeze(dL(3,:,:));
   duLdL = py .* squeeze(dL(2,:,:)) + pz .* squeeze(dL(3,:,:));
   
   % Pick a threshold
   divVthres = 400;
   k = find(duLdL > divVthres & xq<-1.6);
   
   surf(xq,yq,zq,uL,'Linestyle','none'); colorbar;
   caxis([-200 200])
   axis tight equal
   view([-78.3 11.6]); hold on
   scatter3(xq(k),yq(k),zq(k),'r');
   
   dim = [0.2 0.113 0.05 0.02];
   str = sprintf('it=%d, time=%.1fs',...
      filehead.it,filehead.time);
   a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
      'FontWeight','bold','EdgeColor','none');
   
   frame = getframe(gcf);
   
   clf;
   writeVideo(v,frame);

end

v.close
close(4)

%% Movie5

filenameGM = '~/Ganymede/newPIC/run_G8_newPIC/box_FTE_G8_1200s.outs';
filenamePC='~/Ganymede/newPIC/run_G8_newPIC/3d_fluid_35.outs';
%filenamePC='~/Ganymede/newPIC/3d_theta51_test.outs';

v = VideoWriter('E_test3.avi');
v.FrameRate = 10;
v.open

fig = figure(5); 
set(fig,'nextplot','replacechildren');
grid on

for ipict = 1:30
   % GM: clock angle
   [~,data] = read_data(filenameGM,'verbose',false,'npict',ipict);
   
   x = data.file1.x(:,:,:,1);
   y = data.file1.x(:,:,:,2);
   z = data.file1.x(:,:,:,3);
   
   bx = data.file1.w(:,:,:,5);
   by = data.file1.w(:,:,:,6);
   bz = data.file1.w(:,:,:,7);
   
   % From ndgrid to meshgrid format
   x  = permute(x,[2 1 3]);
   y  = permute(y,[2 1 3]);
   z  = permute(z,[2 1 3]);
   bx = permute(bx,[2 1 3]);
   by = permute(by,[2 1 3]);
   bz = permute(bz,[2 1 3]);
   
   % Calculate the inner direction tangential to the boundary surface
   Coef = Coefin;
   bxv= interp3(x, y, z, bx, Coef*xq, Coef*yq, Coef*zq);
   byv= interp3(x, y, z, by, Coef*xq, Coef*yq, Coef*zq);
   bzv= interp3(x, y, z, bz, Coef*xq, Coef*yq, Coef*zq);
   
   % Transform vectors into local coordinate system
   b1in = Inf(size(xq)); b2in = b1in;
   
   % This could potentially be improved!
   for ix=1:size(xq,1)
      for iy=1:size(xq,2)
         b1in(ix,iy)  = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
         b2in(ix,iy)  = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
      end
   end
   
   % Calculate the outer direction tangential to the boundary surface
   Coef = Coefout;
   bxv= interp3(x, y, z, bx, Coef*xq, Coef*yq, Coef*zq);
   byv= interp3(x, y, z, by, Coef*xq, Coef*yq, Coef*zq);
   bzv= interp3(x, y, z, bz, Coef*xq, Coef*yq, Coef*zq);
   
   
   % Transform vectors into local coordinate system
   b1out = Inf(size(xq)); b2out = b1out;
   
   % This could potentially be improved!
   for ix=1:size(xq,1)
      for iy=1:size(xq,2)
         b1out(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
         b2out(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
      end
   end
   
   % Calculate the shear angle
   % equivalent to: acosd( dot(u,v)/(norm(u)*norm(v)) );
   theta = acosd( (b1in.*b1out + b2in.*b2out) ./ ...
      sqrt((b1in.^2+b2in.^2).*(b1out.^2+b2out.^2)) );
   
   
   
   % PC: identifying reconnection site
   [filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);
   
   data = data.file1;
   x = data.x(:,:,:,1);       % [Rg]
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);
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
   
   % Electron pressure gradient + electron intertia term [mV/m]
   UeCrossB = cross([Uex(:) Uey(:) Uez(:)],[Bx(:) By(:) Bz(:)])*1e-3;
   Eprime = [Ex(:) Ey(:) Ez(:)] + UeCrossB;
   
   Eprime = sqrt(sum(Eprime.^2,2));
   Eprime = reshape(Eprime,size(Ex));
   
   % Pick a threshold
   EprimeThres = 60; %[mV/m]
   k = find(Eprime > EprimeThres & z < 0.5 & z > -0.5);
   
   
   % Plotting
   contourf(yq,zq,theta,'Linestyle','none');
   c = colorbar; colormap jet; 
   c.Label.String = 'degree';
   axis equal; caxis manual; caxis([0 180])
   hold on
   scatter(y(k),z(k));
   axis([-1.1 1.1 -0.5 0.5]);
   set(gca,'XDir','reverse')
   
   frame = getframe(gcf);
   
   clf;
   writeVideo(v,frame);

end

v.close
close(5)



%% MovieFinal
% Combine the clock angle and reconnection site identification.

filenameGM = '~/Ganymede/newPIC/run_G8_newPIC/box_FTE_G8_1200s.outs';
filenamePC='~/Ganymede/newPIC/run_G8_newPIC/3d_fluid_35.outs';
%filenamePC='~/Ganymede/newPIC/3d_theta51_test.outs';

v = VideoWriter('E_test2.avi');
v.FrameRate = 10;
v.open

fig = figure(7); 
set(fig,'nextplot','replacechildren');
grid on

for ipict = 1:30
   % GM: clock angle
   [~,data] = read_data(filenameGM,'verbose',false,'npict',ipict);
   
   x = data.file1.x(:,:,:,1);
   y = data.file1.x(:,:,:,2);
   z = data.file1.x(:,:,:,3);
   
   bx = data.file1.w(:,:,:,5);
   by = data.file1.w(:,:,:,6);
   bz = data.file1.w(:,:,:,7);
   
   % From ndgrid to meshgrid format
   x  = permute(x,[2 1 3]);
   y  = permute(y,[2 1 3]);
   z  = permute(z,[2 1 3]);
   bx = permute(bx,[2 1 3]);
   by = permute(by,[2 1 3]);
   bz = permute(bz,[2 1 3]);
   
   % Calculate the inner direction tangential to the boundary surface
   Coef = Coefin;
   bxv= interp3(x, y, z, bx, Coef*xq, Coef*yq, Coef*zq);
   byv= interp3(x, y, z, by, Coef*xq, Coef*yq, Coef*zq);
   bzv= interp3(x, y, z, bz, Coef*xq, Coef*yq, Coef*zq);
   
   % Transform vectors into local coordinate system
   b1in = Inf(size(xq)); b2in = b1in;
   
   % This could potentially be improved!
   for ix=1:size(xq,1)
      for iy=1:size(xq,2)
         b1in(ix,iy)  = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
         b2in(ix,iy)  = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
      end
   end
   
   % Calculate the outer direction tangential to the boundary surface
   Coef = Coefout;
   bxv= interp3(x, y, z, bx, Coef*xq, Coef*yq, Coef*zq);
   byv= interp3(x, y, z, by, Coef*xq, Coef*yq, Coef*zq);
   bzv= interp3(x, y, z, bz, Coef*xq, Coef*yq, Coef*zq);
   
   
   % Transform vectors into local coordinate system
   b1out = Inf(size(xq)); b2out = b1out;
   
   % This could potentially be improved!
   for ix=1:size(xq,1)
      for iy=1:size(xq,2)
         b1out(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
         b2out(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
      end
   end
   
   % Calculate the shear angle
   % equivalent to: acosd( dot(u,v)/(norm(u)*norm(v)) );
   theta = acosd( (b1in.*b1out + b2in.*b2out) ./ ...
      sqrt((b1in.^2+b2in.^2).*(b1out.^2+b2out.^2)) );
   
   
   
   % PC: identifying reconnection site
   [filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);
   
   data = data.file1;
   x = data.x(:,:,:,1);       % [Rg]
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);
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
   
   % Electron pressure gradient + electron intertia term [mV/m]
   UeCrossB = cross([Uex(:) Uey(:) Uez(:)],[Bx(:) By(:) Bz(:)])*1e-3;
   Eprime = [Ex(:) Ey(:) Ez(:)] + UeCrossB;
   
   Eprime = sqrt(sum(Eprime.^2,2));
   Eprime = reshape(Eprime,size(Ex));
   
   % Pick a threshold
   EprimeThres = 60; %[mV/m]
   k = find(Eprime > EprimeThres & z < 0.5 & z > -0.5);
   
   
   % Plotting
   contourf(yq,zq,theta,'Linestyle','none');
   c = colorbar; colormap jet; 
   c.Label.String = 'degree';
   axis equal; caxis manual; caxis([0 180])
   hold on
   scatter(y(k),z(k));
   axis([-1.1 1.1 -0.5 0.5]);
   set(gca,'XDir','reverse')
   
   frame = getframe(gcf);
   
   clf;
   writeVideo(v,frame);

end

v.close
close(7)