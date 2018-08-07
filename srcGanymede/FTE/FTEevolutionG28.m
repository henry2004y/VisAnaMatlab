% Ganymede G28 FTE plots for the paper
% Inherited from the AGU plots.
%
% Plot $B_{normal}$, $j$ and $P_e$ and reconnection efficiency in a
% time-series manner.
%
% Hongyang Zhou, hyzhou@umich.edu  12/04/2017
% Modified on 02/20/2018, 08/07/2018

clear;clc;close all
%% Read data and preprocess
IsGatheredFile = false; % Single or multiple input files
boxdir = '~/Ganymede/MOP2018/runG28_PIC_1200s/GM';

G28 = load('~/Ganymede/MOP2018/ProcessedData/CPCP_G28.mat');
rateG28 = G28.CPCPt ./ G28.Potential_bk;
rateMeanG28  = mean(rateG28);

% box outputs
%filename = '~/Ganymede/newPIC/run_G28_newPIC/box_FTE_G28_1200s.outs';

% Parameters and thresholds
Coef = 1.00; % expansion factor
threshold_pe = 2.1;
threshold_j = 0.52; 

%% Find boundary points from steady state solution
filename3dGM = '~/Ganymede/MOP2018/runG28_PIC_1200s/3d_t=566.out';
s = 0.8; % compact boundary factor [0,1]

[x3bc,y3bc,z3bc] = find_boundary_points( filename3dGM,s );

%% Fit the closed field line boundary with hypersurface

[fitresult,gof] = surface_fit(x3bc,y3bc,z3bc);

%% Generate mesh points from fitted surface
ymin = -1.15; ymax = 1.15; zmin = -0.6; zmax = 0.6;
dy = 1/32; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

% Try to rotate the ndgrid by 15 degrees
% Define rotation angle
theta = 15/180*pi; 
% Define rotation matrix
rot = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 

temp = [yq(:),zq(:)]*rot.' ;
Yrot = reshape(temp(:,1),[numel(yq),1]);
Zrot = reshape(temp(:,2),[numel(yq),1]);

Xrot = fitresult(Yrot,Zrot);

xq = reshape(Xrot,size(yq));
yq = reshape(Yrot,size(yq));
zq = reshape(Zrot,size(yq));

% xq = fitresult(yq,zq);
xq(xq>-1.13) = nan;

% Calculate the normal direction to the fitted surface
[V, W] = differentiate(fitresult, yq, zq);
U = -ones(size(V));

% get the three local directions
% dipole-direction unit vector
unitDipole = [19.26 -16.54 716.8]/sqrt(19.26^2+16.54^2+716.8^2);
% Upstream B unit vector
%unitUpstreamB = [-7 78 -76]/sqrt(7^2+78^2+76^2);

unitL = [0 -sind(15) cosd(15)];

% Initialize local vectors: d1-> M d2->L d3-> N
dL = Inf(3,size(xq,1),size(xq,2));
dM = dL; dN = dL;

% This part could potentially be optimized!
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      dN(:,ix,iy) = [U(ix,iy) V(ix,iy) W(ix,iy)];
      dM(:,ix,iy) = cross(dN(:,ix,iy),unitL);
      dL(:,ix,iy) = cross(dM(:,ix,iy),dN(:,ix,iy));
      
      % Normalization
      dL(:,ix,iy) = dL(:,ix,iy) / norm(dL(:,ix,iy));
      dM(:,ix,iy) = dM(:,ix,iy) / norm(dM(:,ix,iy));
      dN(:,ix,iy) = dN(:,ix,iy) / norm(dN(:,ix,iy));
   end
end

%% Plots from rotated ndgrid in the local coordinate system
if IsGatheredFile
   % single input file case
   [filehead,~,fileinfo] = read_data(filename,'verbose',false);
   npict = fileinfo.npictinfiles; % # of snapshot in the file
   firstpict = 516; dpict = 10; lastpict = firstpict + 2*dpict; 
else
   % multiple input file case
   listing = dir(fullfile(boxdir,'box_var_5*out'));
   [filehead,data] = read_data(...
         fullfile(listing(1).folder,listing(1).name),...
         'verbose',false);
   npict = numel(listing); 
   firstpict = 567; dpict = 15; lastpict = firstpict + 2*dpict; 
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
hfig = figure('position',[10, 10, 800, 660]);
colormap(jet);

iplot = 1;
height = 0.24;

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
   p  = permute(p,[2 1 3]);
   %j  = permute(j,[2 1 3]);
   
   bxv= interp3(x, y, z, bx, Coef*xq, Coef*yq, Coef*zq);
   byv= interp3(x, y, z, by, Coef*xq, Coef*yq, Coef*zq); 
   bzv= interp3(x, y, z, bz, Coef*xq, Coef*yq, Coef*zq);
   jxv= interp3(x, y, z, jx, Coef*xq, Coef*yq, Coef*zq);
   jyv= interp3(x, y, z, jy, Coef*xq, Coef*yq, Coef*zq); 
   jzv= interp3(x, y, z, jz, Coef*xq, Coef*yq, Coef*zq);  
   pev= interp3(x, y, z, pe, Coef*xq, Coef*yq, Coef*zq);
   pv = interp3(x, y, z, p,  Coef*xq, Coef*yq, Coef*zq);
   
   jv = sqrt(jxv.^2 + jyv.^2 + jzv.^2);
      
   % Transform vectors into local coordinate system
   uL = Inf(size(xq)); uM = uL; uN = uL;
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
   
   if iplot == 1
      text(1.3,0,'Bnormal [nT]','rotation',90,...
         'fontsize',14,'FontWeight','bold',...
         'horizontalalignment','center','verticalalignment','bottom');
   elseif iplot == 3
      %h.YTickLabel = [];
      % Get the current axis size
      originalSize = get(gca, 'Position');
      colorbar; 
      set(gca,'Position', originalSize);
   end
   
   axis tight equal   
   h.Position(4) = height;
   %h.XTickLabel = [];   
   h.XDir = 'reverse';
   h.FontSize = 12;
   caxis([-40 40])
   
   % Current density plots
   h = subplot(4,3,iplot+3);
   
   contourf(yq,zq,jv,50,'Linestyle','none'); 
   
   if iplot == 1
      text(1.3,0,'j [\mu A/m^2]','rotation',90,...
         'fontsize',14,'FontWeight','bold',...
         'horizontalalignment','center','verticalalignment','bottom');
   elseif iplot == 2
      xlabel('Y [R_G]','FontWeight','bold','fontsize',14,...
         'horizontalalignment','center','verticalalignment','bottom');
      ylabel('Z [R_G]','FontWeight','bold','fontsize',14,...
         'horizontalalignment','center','verticalalignment','bottom');      
   elseif iplot == 3
      % Get the current axis size
      originalSize = get(gca, 'Position');
      c = colorbar;
      set(gca,'Position', originalSize);
   end
   
   axis tight equal     
   h.Position(4) = height;
   %h.XTickLabel = [];
   h.XDir = 'reverse';
   h.FontSize = 12;
   caxis([0 0.5])
   
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

   axis tight equal   
   h.Position(4) = height;
   h.XDir = 'reverse';
   h.FontSize = 12;   
   caxis([0.1 0.8]);
     
   iplot = iplot + 1;
end

h = subplot(4,3,[10,11,12]);

plot(G28.time,rateG28,'k','linewidth',1.2);
h.Position(4) = height-0.03;
h.FontSize = 14;
vline(firstpict,'b');
vline(firstpict+dpict,'b');
vline(lastpict,'b');
xlabel('time [s]'); 
ylabel('reconnection efficiency')
h.YLim = [0.32,0.52];

% dim = [0.2 0.05 0.05 0.02];
% str = sprintf('it=%d, time=%.1fs, countJ=%d,countP=%d',...
%    filehead.it,filehead.time, FTEcountJ,FTEcountP);
% a = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
%    'FontWeight','bold','EdgeColor','none');


function hhh=vline(x,in1,in2)
% function h=vline(x, linetype, label)
% 
% Draws a vertical line on the current axes at the location specified by 'x'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be 
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = vline(42,'g','The Answer')
%
% returns a handle to a green vertical line on the current axes at x=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% vline also supports vector inputs to draw multiple lines at once.  For example,
%
% vline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
%
% draws three lines with the appropriate labels and colors.
% 
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001

if length(x)>1  % vector input
    for I=1:length(x)
        switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            if ~iscell(in1)
                in1={in1};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            label='';
        case 3
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
            end
        end
        h(I)=vline(x(I),linetype,label);
    end
else
    switch nargin
    case 1
        linetype='r:';
        label='';
    case 2
        linetype=in1;
        label='';
    case 3
        linetype=in1;
        label=in2;
    end
   
    g=ishold(gca);
    hold on

    y=get(gca,'ylim');
    h=plot([x x],y,linetype);
    if ~isempty(label)
        xx=get(gca,'xlim');
        xrange=xx(2)-xx(1);
        xunit=(x-xx(1))/xrange;
        if xunit<0.8
            text(x+0.01*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        else
            text(x-.05*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        end
    end     

    if g==0
    hold off
    end
    set(h,'tag','vline','handlevisibility','off')
end % else

if nargout
    hhh=h;
end

end