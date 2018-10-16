function calc_energy_flux(particle,angle)
%CALC_ENERGY_FLUX Calculate the energetic fluxes at the surface
% First calculate the energetic fluxes along the field line at the original
% locations, and then mapping to the surface.
%
%INPUTS:
% particle: particle positions, velocities and weights, [7,nP]
% angle: pitch angles, [nP,1]

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
m = Parameters.m;

xMin = min(particle(1,:)); xMax = max(particle(1,:));
yMin = min(particle(2,:)); yMax = max(particle(2,:));
zMin = min(particle(3,:)); zMax = max(particle(3,:));

bins = 30;

y = linspace(yMin,yMax,bins);
z = linspace(zMin,zMax,bins);
[Y,Z] = ndgrid(y,z);
dy = (yMax - yMin)/(bins - 1); 
dz = (zMax - zMin)/(bins - 1);
subs = [round((particle(2,:)'-yMin)/dy)+1 ...
        round((particle(3,:)'-zMin)/dz)+1];

vPar = sqrt(sum(particle(4:6,:).^2,1)') .* cosd(angle);
% Volume
Vol = (xMax - xMin)*dy*dz;
     
flux = accumarray(subs,...
   (0.5*m*sum(particle(4:6,:).^2,1)'.*vPar.*particle(7,:)')) ./ Vol;

% figure
% surf(Y',Z',flux')

%% Mapping onto the surface
X = mean(particle(1,:)) * ones(size(Y));
theta1 = Ftheta1(X,Y,Z);
phi1   = Fphi1(X,Y,Z);

% Get B field at the orginal locations
[xF,yF,zF,Bx,By,Bz] = get_field;
Fx = griddedInterpolant(xF,yF,zF,Bx);
Fy = griddedInterpolant(xF,yF,zF,By);
Fz = griddedInterpolant(xF,yF,zF,Bz);

Bx = Fx(X,Y,Z); By = Fy(X,Y,Z); Bz = Fz(X,Y,Z);
B  = sqrt(Bx.^2 + By.^2 + Bz.^2);

[FBxSurf,FBySurf,FBzSurf] = get_Bsurface(false);
BxSurf = FBxSurf(phi1,theta1);
BySurf = FBySurf(phi1,theta1);
BzSurf = FBzSurf(phi1,theta1);
Bsurf = sqrt(BxSurf.^2 + BySurf.^2 + BzSurf.^2);

flux = flux .* Bsurf ./ B; 

% Get the flux normal to the surface
Br = BxSurf .* cosd(phi1) .* cosd(theta1) + ...
     BySurf .* sind(phi1) .* cosd(theta1) + ...
     BzSurf .* sind(theta1);
flux = flux .* abs(Br) ./ Bsurf;
   
% figure
% surf(phi1,theta1,flux)
     
figure
axesm('ortho','origin',[45 45]); 
axis off;
gridm on; 
framem on;
mlabel('equator')
setm(gca,'Origin',[0 180 0])
plabel(120); 
plabel('fontweight','bold')
h1 = surfm(theta1,phi1,flux); %colorbar

end