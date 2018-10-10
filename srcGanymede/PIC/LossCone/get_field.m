function [xF,yF,zF,Bx,By,Bz] = get_field
%GET_PARTICLE Summary of this function goes here
%   Detailed explanation goes here

Dir = '~/Ganymede/MOP2018/runG8_PIC_1200s/EnergeticFlux';

fnameField = '3d_fluid_region0_0_t00000557_n00010710.out';

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

end

