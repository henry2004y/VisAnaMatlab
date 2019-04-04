% Ganymede G8 FTE plots for the paper
% Inherited from the AGU plots.
%
% Plot $B_{normal}$, $j$ and $P_e$ and reconnection efficiency in a
% time-series manner.
%
% Hongyang Zhou, hyzhou@umich.edu  12/04/2017
% Modified on 02/20/2018

clear;clc; close all
%% Read data and preprocess
IsGatheredFile = false; % Single or multiple input files
boxdir = '~/Ganymede/MOP2018/runG8_PIC_1200s/GM';

G8 = load('~/Ganymede/MOP2018/ProcessedData/CPCP_G8.mat');
rateG8 = G8.CPCPt ./ G8.Potential_bk;
rateMeanG8  = mean(rateG8);

% box outputs
%filename = '~/Ganymede/MOP2018/runG8_PIC_1200s/GM/box_FTE_G8_1200s.outs';

% Parameters and thresholds
Coef = 1.00; % expansion factor
threshold_pe = 2.1;
threshold_j = 0.52; 

%% Find boundary points from steady state solution
filename3dGM = '~/Ganymede/MOP2018/runG8_PIC_1200s/3d_t=280.out';
s = 0.8; % compact boundary factor [0,1]

[x3bc,y3bc,z3bc] = find_boundary_points( filename3dGM,s );

%% Fit the closed field line boundary with hypersurface

[fitresult,gof] = surface_fit(x3bc,y3bc,z3bc);

%% Generate mesh points from fitted surface and calculate LMN directions
ymin = -1.2; ymax = 1.2; zmin = -0.6; zmax = 0.6;
dy = 1/32; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

xq = fitresult(yq,zq);

% Calculate the normal direction to the fitted surface
[V, W] = differentiate(fitresult, yq, zq);
U = -ones(size(V));

% get the three local directions
% dipole-direction unit vector
unitDipole = [18 -51.82 716.8]/sqrt(18^2+51.82^2+716.8^2);
% Initialize local vectors: d1-> M d2->L d3-> N
dL = Inf(3,size(xq,1),size(xq,2));
dM = dL; dN = dL;

% This part could potentially be optimized!
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      dN(:,ix,iy) = [U(ix,iy) V(ix,iy) W(ix,iy)];
      dM(:,ix,iy) = cross(dN(:,ix,iy),unitDipole);
      dL(:,ix,iy) = cross(dM(:,ix,iy),dN(:,ix,iy));
      
      % Normalization
      dL(:,ix,iy) = dL(:,ix,iy) / norm(dL(:,ix,iy));
      dM(:,ix,iy) = dM(:,ix,iy) / norm(dM(:,ix,iy));
      dN(:,ix,iy) = dN(:,ix,iy) / norm(dN(:,ix,iy));
   end
end

%% Output variables in the local coordinate system

if IsGatheredFile
   % single input file case
   [filehead,~,fileinfo] = read_data(filename,'verbose',false);
   npict = fileinfo.npictinfiles; % # of snapshot in the file
   firstpict = 750; dpict = 10; lastpict = firstpict + 2*dpict;
else
   % multiple input file case
   listing = dir(fullfile(boxdir,'box_var_5*out'));
   [filehead,data] = read_data(...
         fullfile(listing(1).folder,listing(1).name),...
         'verbose',false);
   npict = numel(listing); 
   firstpict = 300; dpict = 15; lastpict = firstpict + 2*dpict;
end

% Maybe write a function?
n_ = strcmpi('rho',filehead.wnames);
bx_ = strcmpi('bx',filehead.wnames);
by_ = strcmpi('by',filehead.wnames);
bz_ = strcmpi('bz',filehead.wnames);
ux_ = strcmpi('ux',filehead.wnames);
uy_ = strcmpi('uy',filehead.wnames);
uz_ = strcmpi('uz',filehead.wnames);
pe_ = strcmpi('pe',filehead.wnames);
p_ = strcmpi('p',filehead.wnames);
jx_ = strcmpi('jx',filehead.wnames);
jy_ = strcmpi('jy',filehead.wnames);
jz_ = strcmpi('jz',filehead.wnames);


% create new figure with specified size
hfig = figure('position',[10, 10, 800, 520]);
colormap(jet);

iplot = 1;
height = 0.2;

% Loop over snapshots
for ipict = firstpict:dpict:lastpict
   fprintf('ipict=%d\n',ipict)

   if IsGatheredFile
      [filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);
   else
      [filehead,data] = read_data(...
         fullfile(listing(ipict).folder,listing(ipict).name),...
         'verbose',false);
      data = data.file1;
   end
   
   x = data.x(:,:,:,1);
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);
   
   bx = data.w(:,:,:,bx_);
   by = data.w(:,:,:,by_);
   bz = data.w(:,:,:,bz_);
   pe = data.w(:,:,:,pe_);
   p  = data.w(:,:,:,p_);
   jx = data.w(:,:,:,jx_);
   jy = data.w(:,:,:,jy_);
   jz = data.w(:,:,:,jz_);
   %j  = sqrt(jx.^2 + jy.^2 + jz.^2); 
   
   % From ndgrid to meshgrid format
   x  = permute(x,[2 1 3]);
   y  = permute(y,[2 1 3]);
   z  = permute(z,[2 1 3]);
   bx = permute(bx,[2 1 3]);
   by = permute(by,[2 1 3]);
   bz = permute(bz,[2 1 3]);
   jx = permute(jx,[2 1 3]);
   jy = permute(jy,[2 1 3]);
   jz = permute(jz,[2 1 3]); 
   pe = permute(pe,[2 1 3]);
   p  = permute(p ,[2 1 3]);
   %j  = permute(j,[2 1 3]);
      
   bxv= interp3(x, y, z, bx, Coef*xq, Coef*yq, Coef*zq);
   byv= interp3(x, y, z, by, Coef*xq, Coef*yq, Coef*zq); 
   bzv= interp3(x, y, z, bz, Coef*xq, Coef*yq, Coef*zq);
   jxv= interp3(x, y, z, jx, Coef*xq, Coef*yq, Coef*zq);
   jyv= interp3(x, y, z, jy, Coef*xq, Coef*yq, Coef*zq); 
   jzv= interp3(x, y, z, jz, Coef*xq, Coef*yq, Coef*zq);  
   pev= interp3(x, y, z, pe, Coef*xq, Coef*yq, Coef*zq);
   %pv = interp3(x, y, z, p,  Coef*xq, Coef*yq, Coef*zq);  
   
   jv = sqrt(jxv.^2 + jyv.^2 + jzv.^2);
      
   % Transform vectors into local coordinate system
   u1 = Inf(size(xq)); u2 = u1; u3 = u1;
   bL = Inf(size(xq)); bM = bL; bN = bL;
      
   % This could potentially be improved!
   for ix=1:size(xq,1)
      for iy=1:size(xq,2)
         bL(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
         bM(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
         bN(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dN(:,ix,iy);
      end
   end
 
   % Bnormal plots
   h = subplot(4,3,iplot);
   
   contourf(yq,zq,bN,50,'Linestyle','none'); 
   title(sprintf('t=%ds',ipict));
   
   if iplot==1
      text(1.3,0,'Bnormal [nT]','rotation',90,...
         'fontsize',14,'FontWeight','bold',...
         'horizontalalignment','center','verticalalignment','bottom');
   elseif iplot==3
      % Get the current axis size
      originalSize = get(gca, 'Position');
      colorbar; 
      set(gca,'Position', originalSize);
   end
   
   h.Position(4) = height;
   h.XDir = 'reverse';
   h.FontSize = 12;
   axis tight equal
   caxis([-60 60])
   
   % Current density plots
   h = subplot(4,3,iplot+3);
   
   contourf(yq,zq,jv,50,'Linestyle','none'); 
   
   if iplot==1
      text(1.3,0,'j [\mu A/m^2]','rotation',90,...
         'fontsize',14,'FontWeight','bold',...
         'horizontalalignment','center','verticalalignment','bottom');
   elseif iplot==2
      xlabel('Y [R_G]','FontWeight','bold','fontsize',14,...
         'horizontalalignment','center','verticalalignment','bottom');
      ylabel('Z [R_G]','FontWeight','bold','fontsize',14,...
         'horizontalalignment','center','verticalalignment','bottom');      
   elseif iplot==3
      % Get the current axis size
      originalSize = get(gca, 'Position');
      c = colorbar;
      set(gca,'Position', originalSize);
   end
   
   h.Position(4) = height;
   h.XDir = 'reverse';
   h.FontSize = 12;
   axis tight equal
   caxis([0 0.6])
   
   % Electron pressure plots
   h = subplot(4,3,iplot+6);
   
   contourf(yq,zq,pev,50,'Linestyle','none'); 

   if iplot == 1
      text(1.3,0,'P_e [nPa]','rotation',90,...
         'fontsize',14,'FontWeight','bold',...
         'horizontalalignment','center','verticalalignment','bottom');
   elseif iplot == 3
      % Get the current axis size
      originalSize = get(gca, 'Position');
      colorbar;
      set(gca,'Position', originalSize);   
   end   

   h.Position(4) = height;
   h.XDir = 'reverse';
   h.FontSize = 12;   
   axis tight equal
   caxis([0 1.5]);
     
   iplot = iplot + 1;
end


h = subplot(4,3,[10,11,12]);

plot(G8.time,rateG8,'k','linewidth',1.2);
h.Position(4) = height;
h.FontSize = 14;
xline(firstpict,'b');
xline(firstpict+dpict,'b');
xline(lastpict,'b');
xlabel('time [s]'); 
ylabel('reconnection efficiency')
h.YLim = [0.3,0.6];

% dim = [0.2 0.05 0.05 0.02];
% str = sprintf('it=%d, time=%.1fs, countJ=%d,countP=%d',...
%    filehead.it,filehead.time, FTEcountJ,FTEcountP);
% a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
%    'FontWeight','bold','EdgeColor','none');
