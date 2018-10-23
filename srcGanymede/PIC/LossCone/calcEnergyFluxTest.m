function calcEnergyFluxTest(particle)
%CalcEnergyFluxTest Test of flux calculation
%
% Basically I want to check two things:
% 1. magnetic flux conservation B1*A1 = B2*A2;
% 2. flux mapping pattern: is the aurora like brightening due to the
% mapping, or is it because of the energetic particle distribution?
%
%Procedures:
% 1. Set the energy flux density for the original mesh to be uniform;
% 2. Map the flux onto the surface;
% 3. Check the conservation of total energy.
%
%
%
%
%INPUTS:
% particle: particle positions, velocities and weights, [7,nP]

%% Obtain topology for the mesh
Dir = Parameters.Dir;
fnameGM = Parameters.fnameGM;

[filehead,data] = read_data(fullfile(Dir,fnameGM),'verbose',false);
data = data.file1;

theta1_ = strcmpi('theta1',filehead.wnames);
phi1_ = strcmpi('phi1',filehead.wnames);

xGM = data.x(:,:,:,1);       % [Rg]
yGM = data.x(:,:,:,2);       % [Rg]
zGM = data.x(:,:,:,3);       % [Rg]

theta1 = data.w(:,:,:,theta1_);
phi1 = data.w(:,:,:,phi1_);

Ftheta1 = griddedInterpolant(xGM,yGM,zGM,theta1);
Fphi1 = griddedInterpolant(xGM,yGM,zGM,phi1);

clearvars Dir fnameGM theta1_ phi1_ data xGM yGM zGM theta1 phi1

%% Calculate flux at the original locations
xMin = min(particle(1,:)); xMax = max(particle(1,:));
yMin = min(particle(2,:)); yMax = max(particle(2,:));
zMin = min(particle(3,:)); zMax = max(particle(3,:));

bins = Parameters.bins;
dy = (yMax - yMin)/(bins - 1); 
dz = (zMax - zMin)/(bins - 1);

y = linspace(yMin,yMax,bins);
z = linspace(zMin,zMax,bins);
X = mean(particle(1,:)) * ones(bins,bins);
[Y,Z] = ndgrid(y,z);

% Get the B field at cell centers
[xF,yF,zF,Bx,By,Bz] = getField;
% Get B field at grid points
Fx = griddedInterpolant(xF,yF,zF,Bx);
Fy = griddedInterpolant(xF,yF,zF,By);
Fz = griddedInterpolant(xF,yF,zF,Bz);

Bx = Fx(X,Y,Z); By = Fy(X,Y,Z); Bz = Fz(X,Y,Z);

% Get the cell center B field
BxC = Inf(size(Bx)-1); ByC = Inf(size(By)-1); BzC = Inf(size(Bz)-1);
BC = Inf(size(Bx)-1);
for i=1:size(Bx,1)-1
   for j=1:size(Bx,2)-1
      BxC(i,j) = 0.25*(Bx(i,j) + Bx(i+1,j) + Bx(i,j+1) + Bx(i+1,j+1));
      ByC(i,j) = 0.25*(By(i,j) + By(i+1,j) + By(i,j+1) + By(i+1,j+1));
      BzC(i,j) = 0.25*(Bz(i,j) + Bz(i+1,j) + Bz(i,j+1) + Bz(i+1,j+1));
      BC(i,j) = sqrt(BxC(i,j).^2 + ByC(i,j).^2 + BzC(i,j).^2);
   end
end


B = [BxC(:) ByC(:) BzC(:)]; 
dir = repmat([1 0 0],numel(BxC),1);

angle = atan2d(vecnorm(cross(B,dir),2,2),dot(B,dir,2));
angle = reshape(angle,size(BxC));
% Assume flux density is unity along the field line
flux = dy*dz.* cosd(angle);

magFlux1 = flux.*BC;
magFlux1 = sum(magFlux1(:));

%% Mapping onto the surface

% Create cell-centered grid
y = yMin+0.5*dy:dy:yMax-0.5*dy;
z = zMin+0.5*dz:dz:zMax-0.5*dz;
X = mean(particle(1,:)) * ones(bins-1,bins-1);
[Y,Z] = ndgrid(y,z);

theta1C = Ftheta1(X,Y,Z);
phi1C   = Fphi1(X,Y,Z);

% Get B field at the orginal locations
[xF,yF,zF,Bx,By,Bz] = getField;
Fx = griddedInterpolant(xF,yF,zF,Bx);
Fy = griddedInterpolant(xF,yF,zF,By);
Fz = griddedInterpolant(xF,yF,zF,Bz);

Bx = Fx(X,Y,Z); By = Fy(X,Y,Z); Bz = Fz(X,Y,Z);
B  = sqrt(Bx.^2 + By.^2 + Bz.^2);

[FBxSurf,FBySurf,FBzSurf] = getBsurface(false);
BxSurf = FBxSurf(phi1C,theta1C);
BySurf = FBySurf(phi1C,theta1C);
BzSurf = FBzSurf(phi1C,theta1C);
Bsurf = sqrt(BxSurf.^2 + BySurf.^2 + BzSurf.^2);

flux = flux .* Bsurf ./ B; 

figure
contourf(theta1C,phi1C,flux,20)

%%
y = linspace(yMin,yMax,bins);
z = linspace(zMin,zMax,bins);
X = mean(particle(1,:)) * ones(bins,bins);
[Y,Z] = ndgrid(y,z);

theta1 = Ftheta1(X,Y,Z);
phi1   = Fphi1(X,Y,Z);

% surface area at r=1
% Get dtheta*dphi for each cell
addpath('/Users/hyzhou/Documents/MATLAB/surfarea')
[~,areaSurf] = surfarea(phi1/180*pi,theta1/180*pi,zeros(bins,bins));


theta1Center = Inf(bins-1,bins-1);
for i=1:size(theta1,1)-1
   for j=1:size(theta1,2)-1
      theta1Center(i,j) = 0.5*(theta1(i,j) + theta1(i,j+1));
   end
end

% Mapping to the spherical surface
areaSurf = areaSurf .* cosd(theta1C);
% Take B field direction into consideration
Br = BxSurf .* cosd(phi1C) .* cosd(theta1C) + ...
     BySurf .* sind(phi1C) .* cosd(theta1C) + ...
     BzSurf .* sind(theta1C);
flux = flux .* abs(Br) ./ Bsurf;

magFlux2 = areaSurf .* abs(Br);
magFlux2 = sum(magFlux2(:));

figure
axesm('ortho','origin',[45 45]); 
axis off;
gridm on; 
framem on;
mlabel('equator')
setm(gca,'Origin',[0 180 0])
plabel(120); 
plabel('fontweight','bold')
h1 = surfm(theta1,phi1,flux); c = colorbar;
ylabel(c,'[W/m^2]')

% figure
% surf(theta1C,phi1C,flux); c = colorbar;

%%
fprintf('magFlux1 = %f\n',magFlux1);
fprintf('magFlux2 = %f\n',magFlux2);
fprintf('conservation error = %f%%\n',abs(magFlux1-magFlux2)/magFlux1*100);

end

