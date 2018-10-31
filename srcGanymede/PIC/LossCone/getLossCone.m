function [particle,angle,Bsurf,Bx_P,By_P,Bz_P,B_P,theta1_P,phi1_P] = ...
   getLossCone(particle,angle,Bx_P,By_P,Bz_P,B_P,weight)
%GETLOSSCONE Select the particles inside the loss cone
%
%INPUTS
% particle: particle info inside the picked region, [6,nP]
% angle: pitch angles, [1,nP]
% Bx_P,By_P,Bz_P:
% B_P: B field strength at particle locations, [1,nP]
% weight: particle weights, [1,nP] 
%         (note that the raw output contains sign info)
%
%OUTPUTS
% Particle: particle info inside the loss cone, [7,nPLoss]
% vPerp: perpendicular velocities of particles inside the loss cone
% vPar : parallel velocities of particles inside the loss cone
% Bsurf: B field strength at the surface for each particle
% Bx_P,By_P,Bz_P,B_P: B field at particle locations
% phi1_P: footprint phi
% theta1_P: foorprint theta

Dir = Parameters.Dir;
fnameGM = Parameters.fnameGM;

[filehead,data] = read_data(fullfile(Dir,fnameGM),'verbose',false);
data = data.file1;

status_ = strcmpi('status',filehead.wnames);
if strcmp(Parameters.Hemisphere,'north')
   theta_ = strcmpi('theta1',filehead.wnames);
   phi_ = strcmpi('phi1',filehead.wnames);
else
   theta_ = strcmpi('theta2',filehead.wnames);
   phi_ = strcmpi('phi2',filehead.wnames);
end

xGM = data.x(:,:,:,1);       % [Rg]
yGM = data.x(:,:,:,2);       % [Rg]
zGM = data.x(:,:,:,3);       % [Rg]

theta1 = data.w(:,:,:,theta_);
phi1 = data.w(:,:,:,phi_);

Ftheta1 = griddedInterpolant(xGM,yGM,zGM,theta1);
Fphi1 = griddedInterpolant(xGM,yGM,zGM,phi1);

[FBxSurf,FBySurf,FBzSurf] = getBsurface(true);

% Find Bsurface for each particle position that the field connects to
theta1_P = Ftheta1(particle(1,:),particle(2,:),particle(3,:))';
phi1_P = Fphi1(particle(1,:),particle(2,:),particle(3,:))';

BxSurf = FBxSurf(phi1_P,theta1_P);
BySurf = FBySurf(phi1_P,theta1_P);
BzSurf = FBzSurf(phi1_P,theta1_P);

Bsurf = sqrt(BxSurf.^2 + BySurf.^2 + BzSurf.^2);

% Mirror ratios
r_mirror = Bsurf ./ B_P;

% loss cone angle
theta_loss = asind(1./sqrt(r_mirror));

% Particles inside the loss cone
PSelect = angle < theta_loss;
particle = [particle(:,PSelect) ; abs(weight(PSelect)')];
Bsurf  = Bsurf(PSelect);
% BxSurf = BxSurf(PSelect);
% BySurf = BySurf(PSelect);
% BzSurf = BzSurf(PSelect);
Bx_P = Bx_P(PSelect);
By_P = By_P(PSelect);
Bz_P = Bz_P(PSelect);
B_P    = B_P(PSelect);
theta1_P = theta1_P(PSelect);
phi1_P = phi1_P(PSelect);
angle  = angle(PSelect); 

% Plot loss cone angle (weight included)
max_theta_loss = ceil(max(theta_loss));
nbins = max_theta_loss + 1;
[histw] = histwv(theta_loss,abs(weight),0,max_theta_loss,nbins);

figure
bar(linspace(0,max_theta_loss,nbins),histw);
xlabel('Loss cone angle [degree]'); ylabel('#');
set(gca,'FontSize',14);

end
