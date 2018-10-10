% Pitch angle calculation
%
% Hongyang Zhou, hyzhou@umich.edu 10/10/2018

clear; clc; %close all
%% Read data
cLight = 4000;
cAlfven = 253;
Dir = '~/Ganymede/MOP2018/runG8_PIC_1200s/EnergeticFlux';
fnameParticle = 'cut_particles0_region0_1_t00000557_n00010710.out';
fnameField = '3d_fluid_region0_0_t00000557_n00010710.out';

% Particle data
[filehead,data] = read_data(fullfile(Dir,fnameParticle));
data = data.file1;

xP = squeeze(data.x(:,:,:,1));
yP = squeeze(data.x(:,:,:,2));
zP = squeeze(data.x(:,:,:,3));

ux_ = strcmpi('ux',filehead.wnames);
uy_ = strcmpi('uy',filehead.wnames);
uz_ = strcmpi('uz',filehead.wnames);
w_  = strcmpi('weight',filehead.wnames);

ux = data.w(:,:,:,ux_);
uy = data.w(:,:,:,uy_);
uz = data.w(:,:,:,uz_);

% uIndex_ = [find(ux_) find(uy_) find(uz_)];
% uxyz = squeeze(data.w(:,:,:,uIndex_));

weight = squeeze(data.w(:,:,:,w_));

% Field data
[filehead,data] = read_data(fullfile(Dir,fnameField));

data = data.file1;

ne_ = strcmpi('ns0',filehead.wnames);
ni_ = strcmpi('ni0',filehead.wnames);
bx_ = strcmpi('bx',filehead.wnames);
by_ = strcmpi('by',filehead.wnames);
bz_ = strcmpi('bz',filehead.wnames);



xF = data.x(:,:,:,1);       % [Rg]
yF = data.x(:,:,:,2);       % [Rg]
zF = data.x(:,:,:,3);       % [Rg]
ne = data.w(:,:,:,ne_)*1e6;    % [#/m^3]
ni = data.w(:,:,:,ni_)*1e6/14; % [#/m^3] 
Bx = data.w(:,:,:,bx_);      % [nT]
By = data.w(:,:,:,by_);
Bz = data.w(:,:,:,bz_);
% Ex = data.w(:,:,:,6)*1e-3; % [mV/m]
% Ey = data.w(:,:,:,7)*1e-3; % [mV/m]
% Ez = data.w(:,:,:,8)*1e-3; % [mV/m]

% The original data is saved in ndgrid format. For streamline and
% isonormals functions, the input should be in meshgrid format.
% xF  = permute(x,[2 1 3]);

%% Pitch angle calculation

[nP,angle,B_P] = get_pitch_angle(xP,yP,zP,xF,yF,zF,Bx,By,Bz);

Region = [-1.2 -1.125 -2 2 1 2];
ncountmax = 1071080;
particle = Inf(6,ncountmax);

nP = 0;
for iP=1:numel(xP)
   if xP(iP) >= Region(1) && xP(iP) <= Region(2) && ...
      yP(iP) >= Region(3) && yP(iP) <= Region(4) && ...
      zP(iP) >= Region(5) && zP(iP) <= Region(6)
      %particle = [particle; ux(ipar) uy(ipar) uz(ipar)];
      nP = nP + 1;
      particle(:,nP) = [xP(iP) yP(iP) zP(iP) ...
         ux(iP) uy(iP) uz(iP)]';
   end
end

angle = Inf(nP,1);

% Get pitch angle for each particle 
Fx = griddedInterpolant(xF,yF,zF,Bx);
Fy = griddedInterpolant(xF,yF,zF,By);
Fz = griddedInterpolant(xF,yF,zF,Bz);

Bx_P = Inf(nP,1); By_P = Inf(nP,1); Bz_P = Inf(nP,1);
for iP=1:nP
   Bx_P(iP) = Fx(particle(1,iP),particle(2,iP),particle(3,iP));
   By_P(iP) = Fy(particle(1,iP),particle(2,iP),particle(3,iP));
   Bz_P(iP) = Fz(particle(1,iP),particle(2,iP),particle(3,iP));
   
   B = [Bx_P(iP) By_P(iP) Bz_P(iP)]; 
   U = [particle(4,iP) particle(5,iP) particle(6,iP)];
   
   angle(iP) = atan2d(norm(cross(B,U)),dot(B,U));   
end

% B Strength at particle positions
B_P = sqrt(Bx_P.^2 + By_P.^2 + Bz_P.^2);

%% Loss cone 

fnameGM = 'box_var_2_t00000557_n00250489.out';
% Particle data
[filehead,data] = read_data(fullfile(Dir,fnameGM));
data = data.file1;

status_ = strcmpi('status',filehead.wnames);
theta1_ = strcmpi('theta1',filehead.wnames);
phi1_ = strcmpi('phi1',filehead.wnames);
theta2_ = strcmpi('theta2',filehead.wnames);
phi2_ = strcmpi('phi2',filehead.wnames);

xGM = data.x(:,:,:,1);       % [Rg]
yGM = data.x(:,:,:,2);       % [Rg]
zGM = data.x(:,:,:,3);       % [Rg]

theta1 = data.w(:,:,:,theta1_);
phi1 = data.w(:,:,:,phi1_);

Ftheta1 = griddedInterpolant(xGM,yGM,zGM,theta1);
Fphi1 = griddedInterpolant(xGM,yGM,zGM,phi1);


[FBxSurf,FBySurf,FBzSurf] = get_Bsurface(true);

BxSurf = Inf(nP,1); BySurf = Inf(nP,1); BzSurf = Inf(nP,1);
% Find Bsurface for each particle position that the field connects to
for iP=1:nP
   % Find theta1, phi1 for each particle
   theta1_p = Ftheta1(particle(1,iP),particle(2,iP),particle(3,iP));
   phi1_p = Fphi1(particle(1,iP),particle(2,iP),particle(3,iP));
   
   BxSurf(iP) = FBxSurf(phi1_p,theta1_p);
   BySurf(iP) = FBySurf(phi1_p,theta1_p);
   BzSurf(iP) = FBzSurf(phi1_p,theta1_p);
end

return
% Find B at the particle positions



% Mirror ratios
r_mirror = Bsurface / Bpar;

theta_loss = asin(1./sqrt(r_mirror));


return

%% Velocity space plot

% Assume the cut region is small enough s.t. it can be treated as uniform
[dBx,dBy,dBz,limits] = GetMeanField(Dir,fnameParticle,fnameField);

% if not, use interpolation method
%[dBx,dBy,dBz,limits] = GetInterpField(Dir,fnameParticle,fnameField);

% v_perp vs. v_par phase space plot
% approximation: v_par = v_z, v_perp = sqrt(v_x^2 + v_y^2)
dPar = [dBx; dBy; dBz]; % Parallel direction
dPerp1 = cross([0 -1 0]',dPar); % Perpendicular direction in-plane
dPerp2 = cross(dPar,dPerp1); % Perpendicular direction out-of-plane

uPar = uxyz*dPar;
uPerp1 = uxyz*dPerp1;
uPerp2 = uxyz*dPerp2;

figure
%  hist3([uPar,uPerp],'CDataMode','auto','FaceColor','interp',...
%     'Edges',{linspace(-cLight,cLight,50),linspace(0,cLight,50)})
% Xedges = [-Inf linspace(-cLight,cLight,50) Inf];
% Yedges = [-Inf linspace(0,cLight,50) Inf];
% Xedges = [-Inf linspace(-10,10,50) Inf];
% Yedges = [-Inf linspace(0,15,50) Inf];
% h = histogram2(uPar/cAlfven,uPerp/cAlfven,Xedges,Yedges);
h = histogram2(uPar/cAlfven,uPerp1/cAlfven);
% h.XBinLimits = [-5,5];
% h.YBinLimits = [0,5];
h.NumBins = [25 25];
h.Normalization = 'probability';
h.FaceColor = 'flat';

axis equal
xlabel('$u_{\parallel}$','Interpreter','latex')
ylabel('$u_{\perp}$','Interpreter','latex')
colorbar
view(2)
set(gca,'FontSize',16,'LineWidth',1.1)

%% Sample region plot over contour
[filehead,data] = read_data(fullfile(Dir,fnameField),'verbose',false);

% Choose your cut
cut = 'y'; PlaneIndex = 64;
plotrange = [nan nan nan nan];

x = data.file1.x(:,:,:,1);
y = data.file1.x(:,:,:,2);
z = data.file1.x(:,:,:,3);

x = permute(x,[2 1 3]);
y = permute(y,[2 1 3]);
z = permute(z,[2 1 3]);

func = 'Bx'; 
func_ = strcmpi(func,filehead.wnames);
Bx = data.file1.w(:,:,:,func_);
Bx = permute(Bx,[2 1 3]);

func = 'Uxs1'; 
func_ = strcmpi(func,filehead.wnames);
uxs1 = data.file1.w(:,:,:,func_);
uxs1 = permute(uxs1,[2 1 3]);

func = 'Uys1'; 
func_ = strcmpi(func,filehead.wnames);
uys1 = data.file1.w(:,:,:,func_);
uys1 = permute(uys1,[2 1 3]);

func = 'Uzs1'; 
func_ = strcmpi(func,filehead.wnames);
uzs1 = data.file1.w(:,:,:,func_);
uzs1 = permute(uzs1,[2 1 3]);

us1 = sqrt(uxs1.^2 + uys1.^2 + uzs1.^2);

func = 'Bz'; 
func_ = strcmpi(func,filehead.wnames);
Bz = data.file1.w(:,:,:,func_);
Bz = permute(Bz,[2 1 3]);

cut1 = squeeze(x(PlaneIndex,:,:));
cut2 = squeeze(z(PlaneIndex,:,:));
Bx    = squeeze(Bx(PlaneIndex,:,:));
%uxs1    = squeeze(uxs1(PlaneIndex,:,:));
us1    = squeeze(us1(PlaneIndex,:,:));
Bz    = squeeze(Bz(PlaneIndex,:,:));
[~, ~, Bx] = subsurface(cut1, cut2, Bx, plotrange);
%[~, ~, uxs1] = subsurface(cut1, cut2, uxs1, plotrange);
[~, ~, us1] = subsurface(cut1, cut2, us1, plotrange);
[cut1, cut2, Bz] = subsurface(cut1, cut2, Bz, plotrange);

figure
contourf(cut1,cut2,us1,50,'Linestyle','none');
colorbar; axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
title('Ui');
set(gca,'FontSize',14,'LineWidth',1.2)
hold on
% streamline function requires the meshgrid format strictly
s = streamslice(cut1',cut2',Bx',Bz',2,'linear');
for is=1:numel(s)
   s(is).Color = 'w'; % Change streamline color to white
   s(is).LineWidth = 1.5;
end

% Plot the picked cut region
rectangle('Position',[limits(1) limits(5) ...
   limits(2)-limits(1) limits(6)-limits(5)])
