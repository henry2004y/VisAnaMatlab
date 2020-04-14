% Distribution plot from a large PIC domain
%
% Hongyang Zhou, hyzhou@umich.edu 01/30/2020

clear; clc; close all
%% Parameters
% G8
cAlfven = 253;       % average Alfven speed, [km/s]
me = 9.10938356e-31; % electron mass, [kg]
mp = 1.6726219e-27;  % proton mass, [kg]
mi = 14;             % average ion mass [amu]

ParticleType = 'i'; % {'i','e'}
% 1: uy-ux; 2: ux-uz; 3: uy-uz
% 4: v_perpOut-v_perpIn 5: v_perpIn-v_par; 6: v_perpOut-v_par
PlotVType = 6;
%Dir = '/Users/hyzhou/Ganymede/PIC_frontera/Particle/flux_rope';
Dir = '~';

fnameParticleE = 'cut_particles0_region0_1_t00001640_n00020369.out';
fnameParticleI = 'cut_particles1_region0_2_t00001640_n00020369.out';

% fnameParticleE = 'cut_particles0_region0_1_t00001520_n00004093.out';
% fnameParticleI = 'cut_particles1_region0_2_t00001520_n00004093.out';

if ParticleType == 'i'
   fnameParticle = fnameParticleI;
elseif ParticleType == 'e'
   fnameParticle = fnameParticleE;
end

fnameField = strcat('3d_var_region0_0_',fnameParticle(end-22:end));

nBox = 9; % number of box regions

% Define regions
xC = -1.91; % center of boxes
yC = 0.0;
zC = -0.15;
xL = 0.008; % box length in x
yL = 0.2;   % box length in y
zL = 0.03;  % box length in z


%% Classify particles based on locations

Region = cell(nBox,1);
Region{1} = [xC-xL*3/2 xC-xL/2   yC-yL/2 yC+yL/2 zC+zL/2   zC+zL*3/2];
Region{2} = [xC-xL/2   xC+xL/2   yC-yL/2 yC+yL/2 zC+zL/2   zC+zL*3/2];
Region{3} = [xC+xL/2   xC+xL*3/2 yC-yL/2 yC+yL/2 zC+zL/2   zC+zL*3/2];
Region{4} = [xC-xL*3/2 xC-xL/2   yC-yL/2 yC+yL/2 zC-zL/2   zC+zL/2];
Region{5} = [xC-xL/2   xC+xL/2   yC-yL/2 yC+yL/2 zC-zL/2   zC+zL/2];
Region{6} = [xC+xL/2   xC+xL*3/2 yC-yL/2 yC+yL/2 zC-zL/2   zC+zL/2];
Region{7} = [xC-xL*3/2 xC-xL/2   yC-yL/2 yC+yL/2 zC-zL*3/2 zC-zL/2];
Region{8} = [xC-xL/2   xC+xL/2   yC-yL/2 yC+yL/2 zC-zL*3/2 zC-zL/2];
Region{9} = [xC+xL/2   xC+xL*3/2 yC-yL/2 yC+yL/2 zC-zL*3/2 zC-zL/2];


particle = cell(nBox,1);
for i = 1:nBox
   particle{i} = [];
end

[filehead,data,~] = read_data(fullfile(Dir,fnameParticle));
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

for ip = 1:numel(x)
   for iR = 1:nBox
      if x(ip)>=Region{iR}(1) && x(ip)<=Region{iR}(2) && ...
         y(ip)>=Region{iR}(3) && y(ip)<=Region{iR}(4) && ...
         z(ip)>=Region{iR}(5) && z(ip)<=Region{iR}(6)
      
         particle{iR} = [particle{iR}; ux(ip) uy(ip) uz(ip)];
         break
      end
   end
end


%% Velocity distribution plot in 9 regions

figure('Name','Distribution','Position',[1 1 1100 700]);
for iB=1:nBox
   
   if PlotVType <= 3
      uz = particle{iB}(:,3);
      ux = particle{iB}(:,1);
      uy = particle{iB}(:,2);
   else
      [dBx,dBy,dBz] = GetMeanField(Dir,fnameParticle,fnameField,...
         Region{iB});
      
      % v_perp .vs. v_par phase space plot
      % approximation: v_par = v_z, v_perp = sqrt(v_x^2 + v_y^2)
      dPar = [dBx; dBy; dBz]; % Parallel direction
      dPerpI = cross([0 -1 0]',dPar); % Perpendicular direction in-plane
      dPerpO = cross(dPar,dPerpI); % Perpendicular direction out-of-plane
      
      uPar = particle{iB}(:,1:3)*dPar;
      uPerp1 = particle{iB}(:,1:3)*dPerpI;
      uPerp2 = particle{iB}(:,1:3)*dPerpO;
   end
   

   subplot(3,4,iB+ceil(iB/3));
   if PlotVType == 1
      h = histogram2(uy/cAlfven,ux/cAlfven);
   elseif PlotVType == 2
      h = histogram2(ux/cAlfven,uz/cAlfven);
   elseif PlotVType == 3
      h = histogram2(uy/cAlfven,uz/cAlfven);
   elseif PlotVType == 4
      h = histogram2(uPerp2/cAlfven,uPerp1/cAlfven);
   elseif PlotVType == 5
      h = histogram2(uPerp1/cAlfven,uPar/cAlfven);
   elseif PlotVType == 6
      h = histogram2(uPerp2/cAlfven,uPar/cAlfven);
   else
      error('Unknown PlotVType!')
   end
   if ParticleType == 'e'
      h.XBinLimits = [-10,12];
      h.YBinLimits = [-10,10];
      str = 'electron';
   elseif ParticleType == 'i'
      h.XBinLimits = [-3,3];
      h.YBinLimits = [-3,3];
      str = 'ion';
   end
   h.NumBins = [60 60];
   h.Normalization = 'pdf';
   h.FaceColor = 'flat';
   h.EdgeColor = 'none';
   h.ShowEmptyBins='off';   
   axis equal
   
   if PlotVType==1
      xlabel('u_y','FontWeight','bold','FontSize',20)
      ylabel('u_x','FontWeight','bold','FontSize',20)
   elseif PlotVType==2
      xlabel('u_x','FontWeight','bold','FontSize',20)
      ylabel('u_z','FontWeight','bold','FontSize',20)
   elseif PlotVType==3
      xlabel('u_y','FontWeight','bold','FontSize',20)
      ylabel('u_z','FontWeight','bold','FontSize',20)
   elseif PlotVType==4
      xlabel('$u_{\perp 2}$','FontWeight','bold','FontSize',20,'Interpreter','latex')
      ylabel('$u_{\perp 1}$','FontWeight','bold','FontSize',20,'Interpreter','latex')
   elseif PlotVType==5
      xlabel('$u_{\perp 1}$','FontWeight','bold','FontSize',20,'Interpreter','latex')
      ylabel('$u_\parallel$','FontWeight','bold','FontSize',20,'Interpreter','latex')
   elseif PlotVType==6
      xlabel('$u_{\perp 2}$','FontWeight','bold','FontSize',20,'Interpreter','latex')
      ylabel('$u_\parallel$','FontWeight','bold','FontSize',20,'Interpreter','latex') 
   end
   %title(sprintf('PIC region %d',ipict))
   title(sprintf('%d, x[%3.3f,%3.3f], z[%3.3f,%3.3f]',...
      iB,Region{iB}(1:2),Region{iB}(5:6)),'FontSize',18)
   colorbar
   colormap(gca,'hot')
   %caxis([0 100])
   view(2)
   set(gca,'FontSize',10,'LineWidth',1.1, 'ColorScale', 'log')

   text(0.05,0.05,str,'FontSize',14, 'Units','normalized');
end


%% Sample region plot over contour
[filehead,data] = read_data(fullfile(Dir,fnameField),'verbose',false);

% Choose your cut
cut = 'y'; PlaneIndex = 128;
plotrange = [xC-xL*16, xC+xL*16, zC-zL*5, zC+zL*5];

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
Bx   = squeeze(Bx(PlaneIndex,:,:));
Ex   = squeeze(Ex(PlaneIndex,:,:));
Bz   = squeeze(Bz(PlaneIndex,:,:));
[~, ~, Bx] = subsurface(cut1, cut2, Bx, plotrange);
[~, ~, Ex]= subsurface(cut1, cut2, Ex, plotrange);
[cut1, cut2, Bz] = subsurface(cut1, cut2, Bz, plotrange);

%figure('Name','Contour','Position',[1 1 400 700]);
subplot(3,4,[1,5,9])
contourf(cut1,cut2,Ex,50,'Linestyle','none');
colorbar; colormap(gca,brewermap([],'*RdBu')); %colormap(gca,'jet');
caxis([-9e4,9e4])
axis equal; 
xlabel('x [R_G]'); ylabel('z [R_G]');
title('Ex [\mu V/m]');
set(gca,'FontSize',16,'LineWidth',1.2)
hold on
% streamline function requires the meshgrid format strictly
s = streamslice(cut1',cut2',Bx',Bz',1,'linear');
for is=1:numel(s)
   s(is).Color = 'k';
   s(is).LineWidth = 1.3;
end

for iB=1:nBox
   rectangle('Position',[Region{iB}(1) Region{iB}(5) ...
      Region{iB}(2)-Region{iB}(1) ...
      Region{iB}(6)-Region{iB}(5)],'EdgeColor','r','LineWidth',1.5)
end

