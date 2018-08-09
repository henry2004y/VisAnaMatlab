% Distribution Plot together
%
%

clear;clc;close all
%% 
cLight = 4000;
cAlfven = 253;
Dir = '~/Ganymede/MOP2018/runG8_PIC_1200s/Particles/';
fnameParticle = 'cut_particles0_region0_8_t00000040_n00001200.out';
fnameField = '3d_fluid_region0_0_t00000040_n00001200.out';

listing = dir(fullfile(Dir,'cut_particles0*'));
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
%    x = data.x(:,:,:,1); x = permute(x,[2 1 3]);
%    y = data.x(:,:,:,2); y = permute(y,[2 1 3]);
%    z = data.x(:,:,:,3); z = permute(z,[2 1 3]);
   ux = data.w(:,:,:,ux_); ux = permute(ux,[2 1 3]);
   uy = data.w(:,:,:,uy_); uy = permute(uy,[2 1 3]);
   uz = data.w(:,:,:,uz_); uz = permute(uz,[2 1 3]);
   
   % Velocity space plot
   
   [Bx,By,Bz,limits{ipict}] = GetMeanField(Dir,...
      listing(ipict).name,fnameField);
   
   % v_perp .vs. v_par phase space plot
   % approximation: v_par = v_z, v_perp = sqrt(v_x^2 + v_y^2)
   uPar = Inf(numel(ux),1); uPerp = Inf(numel(ux),1);
   for i=1:numel(ux)
      uPar(i) = dot([Bx By Bz],[ux(i) uy(i) uz(i)]);
      uPerp(i)= sqrt(ux(i)*ux(i) + uy(i)*uy(i) + uz(i)*uz(i) - ...
         uPar(i)^2);
   end
   
   figure
   h = histogram2(uPar/cAlfven,uPerp/cAlfven);
   h.XBinLimits = [-5,5];
   h.YBinLimits = [0,5];
   h.NumBins = [25 25];
   h.Normalization = 'probability';
   h.FaceColor = 'flat';
   
   axis equal
   xlabel('$u_{\parallel}$','Interpreter','latex')
   ylabel('$u_{\perp}$','Interpreter','latex')
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
Uzx = data.file1.w(:,:,:,func_);
Uzx = permute(Uzx,[2 1 3]);

func = 'Bz'; 
func_ = strcmpi(func,filehead.wnames);
Bz = data.file1.w(:,:,:,func_);
Bz = permute(Bz,[2 1 3]);

cut1 = squeeze(x(PlaneIndex,:,:));
cut2 = squeeze(z(PlaneIndex,:,:));
Bx    = squeeze(Bx(PlaneIndex,:,:));
Uzx    = squeeze(Uzx(PlaneIndex,:,:));
Bz    = squeeze(Bz(PlaneIndex,:,:));
[~, ~, Bx] = subsurface(cut1, cut2, Bx, plotrange);
[~, ~, Uzx] = subsurface(cut1, cut2, Uzx, plotrange);
[cut1, cut2, Bz] = subsurface(cut1, cut2, Bz, plotrange);

figure(npict+1)
contourf(cut1,cut2,Uzx,50,'Linestyle','none');
colorbar; axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
title('Uzx');
set(gca,'FontSize',14,'LineWidth',1.2)
hold on
% streamline function requires the meshgrid format strictly
s = streamslice(cut1',cut2',Bx',Bz',1,'linear');
for is=1:numel(s)
   s(is).Color = 'w'; % Change streamline color to white
   s(is).LineWidth = 1.5;
end

for ipict=1:npict
   rectangle('Position',[limits{ipict}(1) limits{ipict}(5) ...
      limits{ipict}(2)-limits{ipict}(1) limits{ipict}(6)-limits{ipict}(5)])
end