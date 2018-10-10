function [xP,yP,zP,ux,uy,uz,weight] = get_particle
%GET_PARTICLE Summary of this function goes here
%   Detailed explanation goes here

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

end

