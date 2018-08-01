% Distribution function
%
% Hongyang Zhou, hyzhou@umich.edu 07/27/2018

clear; clc; close all
%% Read data
filename = '~/Documents/research/Ganymede/data/3d_Inside_t=0_electron.out';

[filehead,data,list] = read_data(filename,'verbose',true);

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

%% Phase space distribution
figure
scatter3(x(1:1e3:end),y(1:1e3:end),z(1:1e3:end),'.')
xlabel('x'); ylabel('y'); zlabel('z');
% 
% figure;
% scatter(x(1:1e4:end),uz(1:1e4:end),'.')

figure;
X = [x',uz'];
%hist3(X,'CdataMode','auto')
hist3(X,'CDataMode','auto','FaceColor','interp',...
   'Edges',{linspace(-2,-1.8,50),linspace(-1000,1000,50)})
xlabel('x [R_G]')
ylabel('uz')
colorbar
%view(2)
set(gca,'FontSize',16,'LineWidth',1.1)

figure;
X = [x',ux'];
%hist3(X,'CdataMode','auto')
hist3(X,[50 50],'CDataMode','auto','FaceColor','interp')
xlabel('x [R_G]')
ylabel('ux')
colorbar
%view(2)
set(gca,'FontSize',16,'LineWidth',1.1)

figure;
X = [x',uy'];
%hist3(X,'CdataMode','auto')
hist3(X,[50 50],'CDataMode','auto','FaceColor','interp')
xlabel('x [R_G]')
ylabel('uy')
colorbar
%view(2)
set(gca,'FontSize',16,'LineWidth',1.1)
