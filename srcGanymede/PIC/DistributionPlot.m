% Distribution function
%
% Hongyang Zhou, hyzhou@umich.edu 07/27/2018

clear; clc; %close all
%% Read data
cLight = 4000;
cAlfven = 253;
Dir = '~/Ganymede/MOP2018/runG8_PIC_1200s/Particles/';
fnameParticle = 'cut_particles1_region0_2_t00000040_n00001200.out';
fnameField = '3d_fluid_region0_0_t00000040_n00001200.out';

[filehead,data,list] = read_data(fullfile(Dir,fnameParticle));

x = data.file1.x(:,:,:,1);
y = data.file1.x(:,:,:,2);
z = data.file1.x(:,:,:,3);

x = permute(x,[2 1 3]);
y = permute(y,[2 1 3]);
z = permute(z,[2 1 3]);

ux_ = strcmpi('ux',filehead.wnames);
uy_ = strcmpi('uy',filehead.wnames);
uz_ = strcmpi('uz',filehead.wnames);

ux = data.file1.w(:,:,:,ux_); ux = permute(ux,[2 1 3]);
uy = data.file1.w(:,:,:,uy_); uy = permute(uy,[2 1 3]);
uz = data.file1.w(:,:,:,uz_); uz = permute(uz,[2 1 3]);


%% Phase space distribution, Gabor version
% figure
% scatter3(x(1:1e2:end),y(1:1e2:end),z(1:1e2:end),'.')
% xlabel('x'); ylabel('y'); zlabel('z');
% % 
% % figure;
% % scatter(x(1:1e4:end),uz(1:1e4:end),'.')
% 
% figure;
% X = [x',uz'];
% %hist3(X,'CdataMode','auto')
% hist3(X,'CDataMode','auto','FaceColor','interp',...
%    'Edges',{linspace(-2.2,-2.1,50),linspace(-1000,1000,50)})
% xlabel('x [R_G]')
% ylabel('uz')
% colorbar
% %view(2)
% set(gca,'FontSize',16,'LineWidth',1.1)
% 
% figure;
% X = [x',ux'];
% %hist3(X,'CdataMode','auto')
% hist3(X,[50 50],'CDataMode','auto','FaceColor','interp')
% xlabel('x [R_G]')
% ylabel('ux')
% colorbar
% %view(2)
% set(gca,'FontSize',16,'LineWidth',1.1)
% 
% figure;
% X = [x',uy'];
% %hist3(X,'CdataMode','auto')
% hist3(X,[50 50],'CDataMode','auto','FaceColor','interp')
% xlabel('x [R_G]')
% ylabel('uy')
% colorbar
% %view(2)
% set(gca,'FontSize',16,'LineWidth',1.1)

%% Velocity space plot

[Bx,By,Bz,limits] = GetMeanField(Dir,fnameParticle,fnameField);

% v_perp .vs. v_par phase space plot
% approximation: v_par = v_z, v_perp = sqrt(v_x^2 + v_y^2)
uPar = Inf(numel(ux),1); uPerp = Inf(numel(ux),1);
for i=1:numel(ux)
   uPar(i) = dot([Bx By Bz],[ux(i) uy(i) uz(i)]);
   uPerp(i)= norm([ux(i) uy(i) uz(i)] - uPar(i)*[Bx By Bz]);
end


%X = [uPar,uPerp];
figure
% hist3(X,'CDataMode','auto','FaceColor','interp',...
%    'Edges',{linspace(-cLight,cLight,50),linspace(0,cLight,50)})
% Xedges = [-Inf linspace(-cLight,cLight,50) Inf];
% Yedges = [-Inf linspace(0,cLight,50) Inf];
% Xedges = [-Inf linspace(-10,10,50) Inf];
% Yedges = [-Inf linspace(0,15,50) Inf];
% h = histogram2(uPar/cAlfven,uPerp/cAlfven,Xedges,Yedges);
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



