% Calculate Alfven speed from outputs
%
%

%% Physical Parameters
mu0 = 4*pi*1e-7;    %[H/m]
me  = 9.1094e-31;   %[kg] 
mp  = 1.6726*1e-27; %[kg]
mi  = 14; % average ion mass [amu]
e   = 1.6022e-19; %[C]

%%
PCdir = '~/Ganymede/MOP2018/runG8_PIC_1200s/PC';
PCfile = '';


%%
ne_ = strcmpi('rhos0',filehead.wnames);
ni_ = strcmpi('rhos1',filehead.wnames);
bx_ = strcmpi('bx',filehead.wnames);
by_ = strcmpi('by',filehead.wnames);
bz_ = strcmpi('bz',filehead.wnames);
ne  = data.w(:,:,:,ne_)*1e6;    % [/m^3]
ni  = data.w(:,:,:,ni_)*1e6;    % [/m^3]
bx  = data.w(:,:,:,bx_);      % [nT]
by  = data.w(:,:,:,by_);
bz  = data.w(:,:,:,bz_);

% From ndgrid to meshgrid format
ne = permute(ne,[2 1 3]);
ni = permute(ni,[2 1 3]);

VA = sqrt((bx.^2 + by.^2 + bz.^2) ./ (mu0*(ni*mp+ne*me)))*1e-12;

