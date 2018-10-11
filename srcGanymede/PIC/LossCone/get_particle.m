function [xP,yP,zP,ux,uy,uz,weight] = get_particle(TypeParticle)
%GET_PARTICLE Read particle info from PIC output
%
%INPUT
% TypeParticle: character, {'electron','ion'}
%
%OUTPUTS:
% xP,yP,zP: particle positions
% ux,uy,uz: particle velocities
% weight  : particle weights


Dir = '~/Documents/research/Ganymede/data/EnergeticFlux';
switch TypeParticle
   case 'electron'
      % Electron
%       fnameParticle = 'cut_particles0_region0_1_t00000557_n00010710.out';
      fnameP = Parameters.fnameE;
   case 'ion'
      % Ion
      %fnameParticle = 'cut_particles1_region0_2_t00000557_n00010710.out';
      fnameP = Parameters.fnameI;
   otherwise
      error('unknown particle type!')
end

% Particle data
[filehead,data] = read_data(fullfile(Dir,fnameP),'verbose',false);
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

weight = squeeze(data.w(:,:,:,w_));

end

