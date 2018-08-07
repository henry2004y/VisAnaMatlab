% Cut plot from IPIC3D output
% 
%
% Hongyang Zhou, hyzhou@umich.edu 03/21/2018

clear;clc; close all
%% Example using plot_data
% The cut plot is now incorporated into the plot_data function:
% filename='~/Documents/research/Ganymede/data/3d_fluid.out';
% ipict = 1;
% [filehead,data,list] = read_data(filename,'verbose',true,'npict',ipict);
% plot_data(data.file1,filehead(1),'uzs0','plotrange',[-2 nan -1.2 1.2],...
%    'plotmode','contbar','cut','y','CutPlaneIndex',45);

%%

%filename='~/Ganymede/newPIC/G8_PIC_theta51/3d_fluid_600s.outs';
%ipict = 132;
filename='~/Documents/research/Ganymede/data/3d_fluid.out';
ipict = 1;

% Estimation of Alfven velocity
VA = 253; %[km/s]

[filehead,data,list] = read_data(filename,'verbose',true,'npict',ipict);

% Choose your cut
cut = 'y'; PlaneIndex = 64;
plotrange = [-2.1 -1.7 -1. 1.];
%plotrange = [nan nan nan nan];

x = data.file1.x(:,:,:,1);
y = data.file1.x(:,:,:,2);
z = data.file1.x(:,:,:,3);

x = permute(x,[2 1 3]);
y = permute(y,[2 1 3]);
z = permute(z,[2 1 3]);


%% Magnetic Field

func = 'Bx'; 
func_ = strcmpi(func,filehead.wnames);
Bx = data.file1.w(:,:,:,func_);
Bx = permute(Bx,[2 1 3]);

func = 'By'; 
func_ = strcmpi(func,filehead.wnames);
By = data.file1.w(:,:,:,func_);
By = permute(By,[2 1 3]);

func = 'Bz'; 
func_ = strcmpi(func,filehead.wnames);
Bz = data.file1.w(:,:,:,func_);
Bz = permute(Bz,[2 1 3]);

cut1 = squeeze(x(PlaneIndex,:,:));
cut2 = squeeze(z(PlaneIndex,:,:));
Bx    = squeeze(Bx(PlaneIndex,:,:));
By    = squeeze(By(PlaneIndex,:,:));
Bz    = squeeze(Bz(PlaneIndex,:,:));
[~, ~, Bx] = subsurface(cut1, cut2, Bx, plotrange);
[~, ~, By] = subsurface(cut1, cut2, By, plotrange);
[cut1, cut2, Bz] = subsurface(cut1, cut2, Bz, plotrange);

figure;
contourf(cut1,cut2,By,50,'Linestyle','none');
colorbar; axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
title('By [nT]');
set(gca,'FontSize',14,'LineWidth',1.2)
hold on
% streamline function requires the meshgrid format strictly
s = streamslice(cut1',cut2',Bx',Bz',1,'linear');
for is=1:numel(s)
   s(is).Color = 'w'; % Change streamline color to white
   s(is).LineWidth = 1.5;
end
%quiver(cut1,cut2,Bx,Bz)
%streamline(cut1,cut2,Bx,Bz,-2,2)
%%

func = 'Bx'; 
func_ = strcmpi(func,filehead.wnames);
Bx = data.file1.w(:,:,:,func_);
Bx = permute(Bx,[2 1 3]);

func = 'By'; 
func_ = strcmpi(func,filehead.wnames);
By = data.file1.w(:,:,:,func_);
By = permute(By,[2 1 3]);

func = 'Bz'; 
func_ = strcmpi(func,filehead.wnames);
Bz = data.file1.w(:,:,:,func_);
Bz = permute(Bz,[2 1 3]);
   
figure;
hold on
s = streamslice(x,y,z,Bx,By,Bz,[],0,[]);
for is=1:numel(s)
   s(is).Color = 'k'; % Change streamline color to white
   s(is).LineWidth = 1.5;
end
view(3); axis equal

c = contourslice(x,y,z,By,[],0,[]);
colorbar

% Modify the density of streamlines if needed
% s = streamslice(cut1,cut2,Bx,Bz);
%    
% %    linspace(plotrange(1),plotrange(2),50),...
% %    linspace(plotrange(3),plotrange(4),50));
% 
% for is=1:numel(s)
%    s(is).Color = 'w'; % Change streamline color to white
%    s(is).LineWidth = 1.5;
% end

% hold on
% scatter(cut1(:),cut2(:),'+')
xlabel('x [R_G]'); ylabel('z [R_G]');
%title('By [nT]');
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