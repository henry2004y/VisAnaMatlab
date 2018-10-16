function [xP,yP,zP,ux,uy,uz,weight] = getParticle(TypeParticle)
%GETPARTICLE Read particle info from PIC output
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
      fnameP = Parameters.fnameE;
   case 'ion'
      % Ion
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

ux = data.w(:,:,:,ux_)*1e3; % [m/s]
uy = data.w(:,:,:,uy_)*1e3; % [m/s]
uz = data.w(:,:,:,uz_)*1e3; % [m/s]

weight = squeeze(data.w(:,:,:,w_));

end

