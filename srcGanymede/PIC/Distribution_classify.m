% Distribution plot from a large PIC domain
%
% Hongyang Zhou, hyzhou@umich.edu 08/13/2018

clear; clc; %close all
%% Read data
cAlfven = 253;
Dir = '~/Documents/research/Ganymede/data/DistPlotTest';
fnameParticle = 'cut_particles1_region0_2_t00000710_n00012900.out';
fnameField = '3d_fluid_region0_0_t00000710_n00012900.out';

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
Region{1} = [-2.08 -2.01 -0.05 0.05 0 0.11];
Region{2} = [-1.81 -1.72 -0.08 0.08 0 0.11];
Region{3} = [-1.95 -1.89 -0.08 0.08 -0.14 -0.04];
Region{4} = [-1.86 -1.80 -0.08 0.08 0.36 0.45];

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

%% Distribution plot in four regions

% Velocity space plot
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
   
   figure(1); subplot(2,2,ipict);
   h = histogram2(uPar/cAlfven,uPerp1/cAlfven);
   h.XBinLimits = [-4,4];
   h.YBinLimits = [-4,4];
   h.NumBins = [25 25];
   h.Normalization = 'probability';
   h.FaceColor = 'flat';
   h.EdgeColor = 'none';
   
   axis equal
   xlabel('$u_{\parallel}$','Interpreter','latex')
   ylabel('$u_{\perp1}$','Interpreter','latex')
   colorbar
   view(2)
   set(gca,'FontSize',16,'LineWidth',1.1)

end
