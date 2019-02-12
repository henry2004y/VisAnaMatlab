% Distribution plot from a large PIC domain
%
% Note: The particle weights should be taken into account. Because of
% normalization, the weights, when being transformed to particle mass, is
% acting like real 'weights' for each particle. Thus when plotting the
% phase space distribution, each point in the plane should be multiplied by
% the weight accordingly to get a actual ratio.
%
% I found a function of plotting the weight 2D histogram on FileExchange:
% hist2w
% After some tests, I found that the weighted distribution is very similar
% to the unweighted ones (at least qualitatively), but the former one is
% much slower to run. Therefore I decided to keep the unweighted one for
% now.
%
%
% Hongyang Zhou, hyzhou@umich.edu 08/13/2018

clear; clc; close all
%% Parameters
PlotVType = 2; % 1: v_par vs. v_perp1 2: v_perp1 vs. v_perp2
% G8
Dir = '~/Documents/research/Ganymede/data/DistPlotTest';
fnameParticle = 'cut_particles1_region0_2_t00000710_n00012900.out'; % ion
%fnameParticle = 'cut_particles0_region0_1_t00000710_n00012900.out'; % e
cAlfven = 253; % G8: 253, G28: 450 
me = 9.10938356e-31; % electron mass, [kg]
mp = 1.6726219e-27;  % proton mass, [kg]
mi = 14;             % average ion mass [amu]

TypeParticle = 'ion'; %{'electron','ion'}
PlotVType = 2; % 1: v_par vs. v_perp1 2: v_perp1 vs. v_perp2
Dir = '~/Ganymede/MOP2018/runG8_PIC_1200s/Particles';
fnameParticle = 'cut_particles1_region0_2_t00000710_n00012900.out';
fnameField = '3d_fluid_region0_0_t00000710_n00012900.out';


% G28
%Dir = '~/Documents/research/Ganymede/data/DistPlotTest/G28';
%fnameParticle = 'cut_particles1_region0_2_t00000215_n00005400.out';
%fnameField = '3d_fluid_region0_0_t00000215_n00005400.out';

%% Read data
[filehead,data,list] = read_data(fullfile(Dir,fnameParticle));
data = data.file1;

x = squeeze(data.x(:,:,:,1));
y = squeeze(data.x(:,:,:,2));
z = squeeze(data.x(:,:,:,3));

ux_ = strcmpi('ux',filehead.wnames);
uy_ = strcmpi('uy',filehead.wnames);
uz_ = strcmpi('uz',filehead.wnames);
w_  = strcmpi('weight',filehead.wnames);

ux = data.w(:,:,:,ux_);
uy = data.w(:,:,:,uy_);
uz = data.w(:,:,:,uz_);

% uIndex_ = [find(ux_) find(uy_) find(uz_)];
% uxyz = squeeze(data.w(:,:,:,uIndex_));

switch TypeParticle
   case 'electron'
      % Electron
      mSpecies = 'ms0';
   case 'ion'
      % Ion
      mSpecies = 'ms1';
   otherwise
      error('unknown particle type!')
end

% Obtain the ratio of simulated particle mass to proton mass
ms_ = strcmpi(mSpecies,filehead.variables);
ms_ = circshift(ms_,-numel(filehead.wnames)-filehead.ndim);
weight = squeeze(data.w(:,:,:,w_))*filehead.eqpar(ms_);

% If you want mass in SI unit, multiply weight by No2SiMass from runlog.
% Then to get the number density in SI units, divide by particle mass.
% However, here it is unecessary to do so because in the end it is a
% normalized distribution.

%% Classify particles based on locations
Region = cell(4,1);
Region{1} = [-1.98 -1.95 -0.08 0.08 -0.02 0.09];
Region{2} = [-1.83 -1.80 -0.08 0.08 -0.02 0.09];
Region{3} = [-1.92 -1.89 -0.08 0.08 -0.15 -0.04];
Region{4} = [-1.89 -1.86 -0.08 0.08 0.30 0.41];

% Region{1} = [-1.92 -1.89 -0.08 0.08 -0.01 0.1];
% Region{2} = [-1.89 -1.86 -0.08 0.08 -0.01 0.1];
% Region{3} = [-1.86 -1.83 -0.08 0.08 -0.01 0.1];
% Region{4} = [-1.83 -1.80 -0.08 0.08 -0.02 0.09];

particle = cell(4,1);
particle{1} = [];
particle{2} = [];
particle{3} = [];
particle{4} = [];

for ip = 1:numel(x)
   if x(ip)>=Region{1}(1) && x(ip)<=Region{1}(2) && ...
      y(ip)>=Region{1}(3) && y(ip)<=Region{1}(4) && ...
      z(ip)>=Region{1}(5) && z(ip)<=Region{1}(6)
      particle{1} = [particle{1}; ux(ip) uy(ip) uz(ip) weight(ip)];
   elseif x(ip)>=Region{2}(1) && x(ip)<=Region{2}(2) && ...
      y(ip)>=Region{2}(3) && y(ip)<=Region{2}(4) && ...
      z(ip)>=Region{2}(5) && z(ip)<=Region{2}(6)
      particle{2} = [particle{2}; ux(ip) uy(ip) uz(ip) weight(ip)];
   elseif x(ip)>=Region{3}(1) && x(ip)<=Region{3}(2) && ...
      y(ip)>=Region{3}(3) && y(ip)<=Region{3}(4) && ...
      z(ip)>=Region{3}(5) && z(ip)<=Region{3}(6)
      particle{3} = [particle{3}; ux(ip) uy(ip) uz(ip) weight(ip)];
   elseif x(ip)>=Region{4}(1) && x(ip)<=Region{4}(2) && ...
      y(ip)>=Region{4}(3) && y(ip)<=Region{4}(4) && ...
      z(ip)>=Region{4}(5) && z(ip)<=Region{4}(6)
      particle{4} = [particle{4}; ux(ip) uy(ip) uz(ip) weight(ip)];      
   end
end

%% Velocity distribution plot in four regions

figure('Name','Distribution','Position',[1 1 840 700]);
for ipict=1:4     
   [dBx,dBy,dBz] = GetMeanField(Dir,fnameParticle,fnameField,...
      Region{ipict});
   
   % v_perp .vs. v_par phase space plot
   % approximation: v_par = v_z, v_perp = sqrt(v_x^2 + v_y^2)
   dPar = [dBx; dBy; dBz]; % Parallel direction
   dPerp1 = cross([0 -1 0]',dPar); % Perpendicular direction in-plane
   dPerp2 = cross(dPar,dPerp1); % Perpendicular direction out-of-plane
   
   uPar = particle{ipict}*dPar;
   uPerp1 = particle{ipict}*dPerp1;
   uPerp2 = particle{ipict}*dPerp2;
%    uz = particle{ipict}(:,3);
%    ux = particle{ipict}(:,1);
%    uy = particle{ipict}(:,2);
   
   subplot(2,2,ipict);
   if PlotVType==1
      h = histogram2(uPerp1/cAlfven,uPerp2/cAlfven);
   elseif PlotVType==2
      %h = histogram2(uPerp1/cAlfven,uPerp2/cAlfven);
      h = histogram2(ux/cAlfven,uz/cAlfven);
      %h = histogram2(uy/cAlfven,ux/cAlfven);      
   else
      error('Unknown PlotVType!')
   end
   h.XBinLimits = [-3,3];
   h.YBinLimits = [-3,3];
   h.NumBins = [30 30];
   h.Normalization = 'pdf';
   h.FaceColor = 'flat';
   h.EdgeColor = 'none';
   h.ShowEmptyBins='on';
   
   % This part is used for weighted histogram.
%    weights = particle{ipict}(:,4);
%    c = colormap;
%    binCenterX = linspace(h.XBinLimits(1) + ...
%       0.5*(h.XBinLimits(2)-h.XBinLimits(1))/h.NumBins,...
%       h.XBinLimits(2) - ...
%       0.5*(h.XBinLimits(2)-h.XBinLimits(1))/h.NumBins,...
%       h.NumBins(1));
%    binCenterY = linspace(h.XBinLimits(1) + ...
%       0.5*(h.XBinLimits(2)-h.XBinLimits(1))/h.NumBins,...
%       h.XBinLimits(2) - ...
%       0.5*(h.XBinLimits(2)-h.XBinLimits(1))/h.NumBins,...
%       h.NumBins(1));
%    hist2w([uy/cAlfven ux/cAlfven],weights,binCenterX,binCenterY,c)
   
   axis equal
   
   if PlotVType==1
      xlabel('$u_{\perp 1}$','Interpreter','latex','FontWeight','bold')
      ylabel('$u_{\perp 2}$','Interpreter','latex','FontWeight','bold')
   elseif PlotVType==2
%       xlabel('$u_{\perp1}$','Interpreter','latex')
%       ylabel('$u_{\perp2}$','Interpreter','latex')
      xlabel('u_x','FontWeight','bold','FontSize',20)
      ylabel('u_z','FontWeight','bold','FontSize',20)
   end
   %title(sprintf('PIC region %d',ipict))
   title(sprintf('%d, x[%3.2f,%3.2f], z[%3.2f,%3.2f]',...
      ipict,Region{ipict}(1:2),Region{ipict}(5:6)),'FontSize',18)
   colorbar
   view(2)
   set(gca,'FontSize',18,'LineWidth',1.1)
end


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

func = 'Ex'; 
func_ = strcmpi(func,filehead.wnames);
Ex = data.file1.w(:,:,:,func_);
Ex = permute(Ex,[2 1 3]);

func = 'Bz'; 
func_ = strcmpi(func,filehead.wnames);
Bz = data.file1.w(:,:,:,func_);
Bz = permute(Bz,[2 1 3]);

cut1 = squeeze(x(PlaneIndex,:,:));
cut2 = squeeze(z(PlaneIndex,:,:));
Bx    = squeeze(Bx(PlaneIndex,:,:));
Ex   = squeeze(Ex(PlaneIndex,:,:));
Bz    = squeeze(Bz(PlaneIndex,:,:));
[~, ~, Bx] = subsurface(cut1, cut2, Bx, plotrange);
[~, ~, Ex]= subsurface(cut1, cut2, Ex, plotrange);
[cut1, cut2, Bz] = subsurface(cut1, cut2, Bz, plotrange);

figure('Name','Contour','Position',[1 1 400 700]);
contourf(cut1,cut2,Ex,50,'Linestyle','none');
colorbar; axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
title('Ex [\mu V/m]');
set(gca,'FontSize',16,'LineWidth',1.2)
hold on
% streamline function requires the meshgrid format strictly
s = streamslice(cut1',cut2',Bx',Bz',5,'linear');
for is=1:numel(s)
   s(is).Color = 'w'; % Change streamline color to white
   s(is).LineWidth = 1.5;
end

for ipict=1:4
   rectangle('Position',[Region{ipict}(1) Region{ipict}(5) ...
      Region{ipict}(2)-Region{ipict}(1) ...
      Region{ipict}(6)-Region{ipict}(5)],'EdgeColor','r','LineWidth',1.5)
end

