% Cut plots from PIC 3D outputs for Ganymede paper
%
% Hall magnetic field: By
% Electron outflow velocity: uzs1
% Electron drift velocity: uys1
% Ion velocity?
%
% The cut plot is now incorporated into the plot_data function:
% 
% Example:
% filename = '../newPIC/highResoTest/3d_64_test.out';
% ipict = 1;
% [filehead,data,list] = read_data(filename,'verbose',true,'npict',ipict);
% plot_data(data.file1,filehead(1),'uzs0','plotrange',[-2 nan -1.2 1.2],...
%    'plotmode','contbar','cut','y','CutPlaneIndex',45);
%
%
% Hongyang Zhou, hyzhou@umich.edu 03/21/2018

clear;clc; %close all
%%

%filename='~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600s.outs';
%ipict = 132;
%filename='~/Ganymede/MOP2018/runG8_PIC_1200s/PC/3d_var_region0_0_t00000520_n00009629.out';
filename='~/Ganymede/MOP2018/runG8_PIC_1200s/PC/3d_var_region0_0_t00000450_n00008729.out';
ipict = 1;

% Estimation of Alfven velocity
VA = 253; %[km/s]

[filehead,data,list] = read_data(filename,'verbose',true,'npict',ipict);

% Choose your cut
cut = 'y'; PlaneIndex = 64;
plotrange = [-2.2 -1.5 -1.2 1.2];

x = data.file1.x(:,:,:,1);
y = data.file1.x(:,:,:,2);
z = data.file1.x(:,:,:,3);

x = permute(x,[2 1 3]);
y = permute(y,[2 1 3]);
z = permute(z,[2 1 3]);


%% Hall Magnetic Field

func = 'By'; 
func_ = strcmpi(func,filehead.wnames);

w = data.file1.w(:,:,:,func_);
w = permute(w,[2 1 3]);


switch cut
   case 'x'
      cut1 = squeeze(y(:,PlaneIndex,:));
      cut2 = squeeze(z(:,PlaneIndex,:));
      w    = squeeze(w(:,PlaneIndex,:));      
   case 'y'
      cut1 = squeeze(x(PlaneIndex,:,:));
      cut2 = squeeze(z(PlaneIndex,:,:));
      w    = squeeze(w(PlaneIndex,:,:));
      [cut1, cut2, w] = subsurface(cut1, cut2, w, plotrange);
   case 'z'
      cut1 = squeeze(x(:,:,PlaneIndex));
      cut2 = squeeze(y(:,:,PlaneIndex));
      w    = squeeze(w(:,:,PlaneIndex));      
end


figure(1);
subplot(231)
contourf(cut1,cut2,w,50,'Linestyle','none');
colorbar; axis equal; 
% hold on
% scatter(cut1(:),cut2(:),'+')
xlabel('x [R_G]'); ylabel('z [R_G]');
title('(a) By [nT]');
set(gca,'FontSize',14,'LineWidth',1.2)

%% Electron Velocity uzs0

func = 'uzs0'; 
func_ = strcmpi(func,filehead.wnames);

w = data.file1.w(:,:,:,func_);
w = permute(w,[2 1 3]);

switch cut
   case 'x'
      cut1 = squeeze(y(:,PlaneIndex,:));
      cut2 = squeeze(z(:,PlaneIndex,:));
      w    = squeeze(w(:,PlaneIndex,:));      
   case 'y'
      cut1 = squeeze(x(PlaneIndex,:,:));
      cut2 = squeeze(z(PlaneIndex,:,:));
      w    = squeeze(w(PlaneIndex,:,:));
      [cut1, cut2, w] = subsurface(cut1, cut2, w, plotrange);
   case 'z'
      cut1 = squeeze(x(:,:,PlaneIndex));
      cut2 = squeeze(y(:,:,PlaneIndex));
      w    = squeeze(w(:,:,PlaneIndex));      
end

figure(1);
subplot(235)
contourf(cut1,cut2,w./VA,50,'Linestyle','none');
colorbar; axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
title('(e) U_{ez}');
set(gca,'FontSize',14,'LineWidth',1.2)


%% Electron Velocity uys0

func = 'uys0'; 
func_ = strcmpi(func,filehead.wnames);

w = data.file1.w(:,:,:,func_);
w = permute(w,[2 1 3]);

switch cut
   case 'x'
      cut1 = squeeze(y(:,PlaneIndex,:));
      cut2 = squeeze(z(:,PlaneIndex,:));
      w    = squeeze(w(:,PlaneIndex,:));      
   case 'y'
      cut1 = squeeze(x(PlaneIndex,:,:));
      cut2 = squeeze(z(PlaneIndex,:,:));
      w    = squeeze(w(PlaneIndex,:,:));
      [cut1, cut2, w] = subsurface(cut1, cut2, w, plotrange);
   case 'z'
      cut1 = squeeze(x(:,:,PlaneIndex));
      cut2 = squeeze(y(:,:,PlaneIndex));
      w    = squeeze(w(:,:,PlaneIndex));      
end

figure(1);
subplot(234)
contourf(cut1,cut2,w./VA,50,'Linestyle','none');
colorbar; axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
title('(d) U_{ey}');
set(gca,'FontSize',14,'LineWidth',1.2)

%% Ion Velocity uzs1

func = 'uzs1'; 
func_ = strcmpi(func,filehead.wnames);

w = data.file1.w(:,:,:,func_);
w = permute(w,[2 1 3]);

switch cut
   case 'x'
      cut1 = squeeze(y(:,PlaneIndex,:));
      cut2 = squeeze(z(:,PlaneIndex,:));
      w    = squeeze(w(:,PlaneIndex,:));      
   case 'y'
      cut1 = squeeze(x(PlaneIndex,:,:));
      cut2 = squeeze(z(PlaneIndex,:,:));
      w    = squeeze(w(PlaneIndex,:,:));
      [cut1, cut2, w] = subsurface(cut1, cut2, w, plotrange);
   case 'z'
      cut1 = squeeze(x(:,:,PlaneIndex));
      cut2 = squeeze(y(:,:,PlaneIndex));
      w    = squeeze(w(:,:,PlaneIndex));      
end

figure(1);
subplot(236)
contourf(cut1,cut2,w./VA,50,'Linestyle','none');
colorbar; axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
title('(f) U_{iz}');
set(gca,'FontSize',14,'LineWidth',1.2)

%% Electron Pressure

func = 'Ps0'; 
func_ = strcmpi(func,filehead.wnames);

w = data.file1.w(:,:,:,func_);
w = permute(w,[2 1 3]);

switch cut
   case 'x'
      cut1 = squeeze(y(:,PlaneIndex,:));
      cut2 = squeeze(z(:,PlaneIndex,:));
      w    = squeeze(w(:,PlaneIndex,:));      
   case 'y'
      cut1 = squeeze(x(PlaneIndex,:,:));
      cut2 = squeeze(z(PlaneIndex,:,:));
      w    = squeeze(w(PlaneIndex,:,:));
      [cut1, cut2, w] = subsurface(cut1, cut2, w, plotrange);
   case 'z'
      cut1 = squeeze(x(:,:,PlaneIndex));
      cut2 = squeeze(y(:,:,PlaneIndex));
      w    = squeeze(w(:,:,PlaneIndex));      
end

figure(1)
subplot(232)
contourf(cut1,cut2,w,50,'Linestyle','none');
colorbar; axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
title('(b) P_{e} [nPa]');
set(gca,'FontSize',14,'LineWidth',1.2)

%% Ion Pressure

func = 'Ps1'; 
func_ = strcmpi(func,filehead.wnames);

w = data.file1.w(:,:,:,func_);
w = permute(w,[2 1 3]);

switch cut
   case 'x'
      cut1 = squeeze(y(:,PlaneIndex,:));
      cut2 = squeeze(z(:,PlaneIndex,:));
      w    = squeeze(w(:,PlaneIndex,:));      
   case 'y'
      cut1 = squeeze(x(PlaneIndex,:,:));
      cut2 = squeeze(z(PlaneIndex,:,:));
      w    = squeeze(w(PlaneIndex,:,:));
      [cut1, cut2, w] = subsurface(cut1, cut2, w, plotrange);
   case 'z'
      cut1 = squeeze(x(:,:,PlaneIndex));
      cut2 = squeeze(y(:,:,PlaneIndex));
      w    = squeeze(w(:,:,PlaneIndex));      
end

figure(1)
subplot(233)
contourf(cut1,cut2,w,50,'Linestyle','none');
colorbar; axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
title('(c) P_{i} [nPa]');
set(gca,'FontSize',14,'LineWidth',1.2)


%%

function [newx, newy, newdata] = subsurface(varargin)
%SUBSURFACE Extract subset of surface dataset.
%  This is a simplified version of subvolume.

x      = varargin{1};
y      = varargin{2};
data   = varargin{3};
limits = varargin{4};

if numel(limits)~=4
  error('Reduction must be [xmin xmax ymin ymax]');
end

if limits(1) > limits(2)
  error(message('MATLAB:subvolume:InvalidReductionXRange'));
end
if limits(3) > limits(4)
  error(message('MATLAB:subvolume:InvalidReductionYRange'));
end

sz = size(data);

hx = x(:,1);
hy = y(1,:);

if isnan(limits(1)),  limits(1) = min(hx); end
if isnan(limits(3)),  limits(3) = min(hy); end
if isnan(limits(2)),  limits(2) = max(hx); end
if isnan(limits(4)),  limits(4) = max(hy); end

xind = find(limits(1)<=hx & hx<=limits(2));
yind = find(limits(3)<=hy & hy<=limits(4));

newdata = subdata(data, xind, yind, sz);

newx = x(xind, yind);
newy = y(xind, yind);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newdata = subdata(data, xind, yind, sz)
newdata = data(xind, yind);
newsz = size(newdata);
if length(sz)>2
  newdata = reshape(newdata, [newsz(1:3) sz(4:end)]);
end

end

end