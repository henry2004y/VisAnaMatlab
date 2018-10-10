% Distribution Plot together
%
%

clear;clc;close all
%% 
cLight = 4000;
cAlfven = 253;
Dir = '~/Documents/research/Ganymede/data/DistPlotTest';
fnameParticle = 'cut_particles0_region0_8_t00000040_n00001200.out';
fnameField = '3d_fluid_region0_0_t00000040_n00001200.out';

listing = dir(fullfile(Dir,'cut_particles0*t00000040*'));
npict = numel(listing);

[filehead,~] = read_data(...
         fullfile(listing(1).folder,listing(1).name),'verbose',false);

ux_ = strcmpi('ux',filehead.wnames);
uy_ = strcmpi('uy',filehead.wnames);
uz_ = strcmpi('uz',filehead.wnames);
      
limits = cell(4,1);

for ipict=1:npict
   [~,data] = read_data(...
         fullfile(listing(ipict).folder,listing(ipict).name),...
         'verbose',false);

   data = data.file1;
   uxyz = squeeze(data.w(:,:,:,1:3));
   
   % Velocity space plot
   
   [dBx,dBy,dBz,limits{ipict}] = GetMeanField(Dir,...
      listing(ipict).name,fnameField);
   
   % v_perp .vs. v_par phase space plot
   % approximation: v_par = v_z, v_perp = sqrt(v_x^2 + v_y^2)
   dPar = [dBx; dBy; dBz]; % Parallel direction
   dPerp1 = cross([0 -1 0]',dPar); % Perpendicular direction in-plane
   dPerp2 = cross(dPar,dPerp1); % Perpendicular direction out-of-plane
   
   uPar = uxyz*dPar;
   uPerp1 = uxyz*dPerp1;
   uPerp2 = uxyz*dPerp2;
   
   figure; subplot(2,2,ipict);
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

func = 'Uzs0'; 
func_ = strcmpi(func,filehead.wnames);
Uez = data.file1.w(:,:,:,func_);
Uez = permute(Uez,[2 1 3]);

func = 'Bz'; 
func_ = strcmpi(func,filehead.wnames);
Bz = data.file1.w(:,:,:,func_);
Bz = permute(Bz,[2 1 3]);

cut1 = squeeze(x(PlaneIndex,:,:));
cut2 = squeeze(z(PlaneIndex,:,:));
Bx    = squeeze(Bx(PlaneIndex,:,:));
Uez    = squeeze(Uez(PlaneIndex,:,:));
Bz    = squeeze(Bz(PlaneIndex,:,:));
[~, ~, Bx] = subsurface(cut1, cut2, Bx, plotrange);
[~, ~, Uez] = subsurface(cut1, cut2, Uez, plotrange);
[cut1, cut2, Bz] = subsurface(cut1, cut2, Bz, plotrange);

figure(npict+1)
contourf(cut1,cut2,Uez,50,'Linestyle','none');
colorbar; axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
title('Uez');
set(gca,'FontSize',14,'LineWidth',1.2)
hold on
% streamline function requires the meshgrid format strictly
s = streamslice(cut1',cut2',Bx',Bz',2,'linear');
for is=1:numel(s)
   s(is).Color = 'w'; % Change streamline color to white
   s(is).LineWidth = 1.5;
end

for ipict=1:npict
   rectangle('Position',[limits{ipict}(1) limits{ipict}(5) ...
      limits{ipict}(2)-limits{ipict}(1) limits{ipict}(6)-limits{ipict}(5)])
end