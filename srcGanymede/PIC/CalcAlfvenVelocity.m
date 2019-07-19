% Calculate Alfven speed from outputs
%
%

clear; clc; close all
%% Physical Parameters
mu0 = 4*pi*1e-7;    %[H/m]
me  = 9.1094e-31;   %[kg] 
mp  = 1.6726*1e-27; %[kg]
mi  = 14; % average ion mass [amu]
e   = 1.6022e-19; %[C]

%%
PCdir = '~/Ganymede/MOP2018/runG8_PIC_1200s/PC';
PCfile = '3d_var_region0_0_t00000926_n00016979.out';

[filehead,data] = read_data(fullfile(PCdir,PCfile),'verbose',false);
data = data.file1;

%%
ne_ = strcmpi('rhos0',filehead.wnames);
ni_ = strcmpi('rhos1',filehead.wnames);
bx_ = strcmpi('bx',filehead.wnames);
by_ = strcmpi('by',filehead.wnames);
bz_ = strcmpi('bz',filehead.wnames);
uix_ = strcmpi('uxs1',filehead.wnames);
uiy_ = strcmpi('uys1',filehead.wnames);
uiz_ = strcmpi('uzs1',filehead.wnames);
x = data.x(:,:,:,1);
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);
ne  = data.w(:,:,:,ne_);    % [cc]
ni  = data.w(:,:,:,ni_);    % [cc]
bx  = data.w(:,:,:,bx_);      % [nT]
by  = data.w(:,:,:,by_);
bz  = data.w(:,:,:,bz_);
uix  = data.w(:,:,:,uix_);
uiy  = data.w(:,:,:,uiy_);
uiz  = data.w(:,:,:,uiz_);

% From ndgrid to meshgrid format
x  = permute(x,[2 1 3]);
y  = permute(y,[2 1 3]);
z  = permute(z,[2 1 3]);
ne = permute(ne,[2 1 3]);
ni = permute(ni,[2 1 3]);
bx = permute(bx,[2 1 3]);
by = permute(by,[2 1 3]);
bz = permute(bz,[2 1 3]);
uix = permute(uix,[2 1 3]);
uiy = permute(uiy,[2 1 3]);
uiz = permute(uiz,[2 1 3]);

VA = sqrt((bx.^2 + by.^2 + bz.^2) ./ (mu0*(ni*mp+ne*me)))*1e-15;

x = squeeze(x(64,:,:));
y = squeeze(y(64,:,:));
z = squeeze(z(64,:,:));
VA = squeeze(VA(64,:,:));
ne = squeeze(ne(64,:,:));
ni = squeeze(ni(64,:,:));

figure
contourf(x,z,VA); colorbar
axis equal
% hold on
% xLineCenter = -1.906;
% plot([-0.6+xLineCenter 0.6+xLineCenter],[0 0],'--k')

figure
contourf(x,z,ni+ne); colorbar
axis equal 

ui = sqrt(uix.^2 + uiy.^2 + uiz.^2);


% c = contourslice(x,y,z,VA,[],[],0,100); axis equal; colorbar
% xlabel('x'); ylabel('y'); zlabel('z');
% 
% figure;
% c = contourslice(x,y,z,VA,[],0,[],100); axis equal; colorbar
% xlabel('x'); ylabel('y'); zlabel('z');