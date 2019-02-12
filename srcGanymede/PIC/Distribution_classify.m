% Distribution plot from a large PIC domain
%
% Hongyang Zhou, hyzhou@umich.edu 08/13/2018

clear; clc; close all
%% Parameters
%cAlfven = 450; % G28
%cAlfven = 253; % G8


PlotVType = 2; % 1: v_par vs. v_perp1 2: v_perp1 vs. v_perp2
% G8
Dir = '~/Documents/research/Ganymede/data/DistPlotTest';
fnameParticle = 'cut_particles1_region0_2_t00000710_n00012900.out'; % ion
%fnameParticle = 'cut_particles0_region0_1_t00000710_n00012900.out'; % e
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

ux = data.w(:,:,:,ux_);
uy = data.w(:,:,:,uy_);
uz = data.w(:,:,:,uz_);

% uIndex_ = [find(ux_) find(uy_) find(uz_)];
% uxyz = squeeze(data.w(:,:,:,uIndex_));

%% Classify particles based on locations
Region = cell(4,1);
Region{1} = [-1.98 -1.95 -0.08 0.08 -0.02 0.09];
Region{2} = [-1.83 -1.80 -0.08 0.08 -0.02 0.09];
Region{3} = [-1.92 -1.89 -0.08 0.08 -0.15 -0.04];
Region{4} = [-1.89 -1.86 -0.08 0.08 0.30 0.41];


particle = cell(4,1);
particle{1} = [];
particle{2} = [];
particle{3} = [];
particle{4} = [];

for ipar = 1:numel(x)
   if x(ipar)>=Region{1}(1) && x(ipar)<=Region{1}(2) && ...
      y(ipar)>=Region{1}(3) && y(ipar)<=Region{1}(4) && ...
      z(ipar)>=Region{1}(5) && z(ipar)<=Region{1}(6)
      particle{1} = [particle{1}; ux(ipar) uy(ipar) uz(ipar)];
   elseif x(ipar)>=Region{2}(1) && x(ipar)<=Region{2}(2) && ...
      y(ipar)>=Region{2}(3) && y(ipar)<=Region{2}(4) && ...
      z(ipar)>=Region{2}(5) && z(ipar)<=Region{2}(6)
      particle{2} = [particle{2}; ux(ipar) uy(ipar) uz(ipar)];
   elseif x(ipar)>=Region{3}(1) && x(ipar)<=Region{3}(2) && ...
      y(ipar)>=Region{3}(3) && y(ipar)<=Region{3}(4) && ...
      z(ipar)>=Region{3}(5) && z(ipar)<=Region{3}(6)
      particle{3} = [particle{3}; ux(ipar) uy(ipar) uz(ipar)];
   elseif x(ipar)>=Region{4}(1) && x(ipar)<=Region{4}(2) && ...
      y(ipar)>=Region{4}(3) && y(ipar)<=Region{4}(4) && ...
      z(ipar)>=Region{4}(5) && z(ipar)<=Region{4}(6)
      particle{4} = [particle{4}; ux(ipar) uy(ipar) uz(ipar)];      
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
      h = histogram2(uPerp1/cAlfven,uPerp2/cAlfven);
      %h = histogram2(ux/cAlfven,uz/cAlfven);
      %h = histogram2(uy/cAlfven,ux/cAlfven);
      %h = histogram2(ux/cAlfven,(uz.*cosd(45)+uy.*cosd(45))./cAlfven);
   else
      error('Unknown PlotVType!')
   end
   h.XBinLimits = [-3,3];
   h.YBinLimits = [-3,3];
   h.NumBins = [30 30];
   h.Normalization = 'probability';
   h.FaceColor = 'flat';
   h.EdgeColor = 'none';
   
   axis equal
   if PlotVType==1
      xlabel('$u_{\perp 1}$','Interpreter','latex','FontWeight','bold')
      ylabel('$u_{\perp 2}$','Interpreter','latex','FontWeight','bold')
   elseif PlotVType==2
%       xlabel('$u_{\perp1}$','Interpreter','latex')
%       ylabel('$u_{\perp2}$','Interpreter','latex')
      xlabel('$u_x$','Interpreter','latex','FontWeight','bold')
      ylabel('$u_z$','Interpreter','latex','FontWeight','bold')
   end
   title(sprintf('PIC region %d',ipict))
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
Ey = data.file1.w(:,:,:,func_);
Ey = permute(Ey,[2 1 3]);

func = 'Bz'; 
func_ = strcmpi(func,filehead.wnames);
Bz = data.file1.w(:,:,:,func_);
Bz = permute(Bz,[2 1 3]);

cut1 = squeeze(x(PlaneIndex,:,:));
cut2 = squeeze(z(PlaneIndex,:,:));
Bx    = squeeze(Bx(PlaneIndex,:,:));
Ey    = squeeze(Ey(PlaneIndex,:,:));
Bz    = squeeze(Bz(PlaneIndex,:,:));
[~, ~, Bx] = subsurface(cut1, cut2, Bx, plotrange);
[~, ~, Ey] = subsurface(cut1, cut2, Ey, plotrange);
[cut1, cut2, Bz] = subsurface(cut1, cut2, Bz, plotrange);

figure(2)
contourf(cut1,cut2,Ey,50,'Linestyle','none');
colorbar; axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
title('Ex [\mu V/m]');
%title('U_{ey} [km/s]')
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

