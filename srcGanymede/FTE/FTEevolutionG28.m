% Ganymede G28 FTE plots for the paper
% Inherited from the AGU plots.
%
% Plot $B_{normal}$, $j$ and $P_e$ and reconnection efficiency in a
% time-series manner.
%
% Hongyang Zhou, hyzhou@umich.edu  12/04/2017
% Modified on 02/20/2018

clear;clc
%% Read data and preprocess
G28 = load('~/Ganymede/newPIC/CPCP_G28.mat');
rateG28 = G28.CPCPt ./ G28.Potential_bk;
rateMeanG28  = mean(rateG28);

%% Find boundary points from steady state solution
filename = '~/Ganymede/newPIC/run_G28_newPIC/3d_G28_steady.out';
s = 0.5; % compact boundary factor [0,1]

[x3bc,y3bc,z3bc] = find_boundary_points( filename,s );

%% Fit the closed field line boundary

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

%% Generate mesh points from fitted surface
ymin = -1.1+5/30; ymax = 1.1-5/30; zmin = -0.5; zmax = 0.75-3/30;
dy = 1/30; dz = dy;
[yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

% Try to rotate the ndgrid by 11 degrees
% Define rotation angle
theta = 11/180*pi; 
% Define rotation matrix
rot = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 

temp = [yq(:),zq(:)]*rot.' ;
sz = numel(yq);
Yq = reshape(temp(:,1),[sz,1]);
Zq = reshape(temp(:,2),[sz,1]);
Xq = fitresult(Yq,Zq);

Xq = reshape(Xq,size(yq));
Yq = reshape(Yq,size(yq));
Zq = reshape(Zq,size(yq));

%% get the three local directions
z = [0 0 1]; % z-direction unit vector
% Initialize local vectors
dL = Inf(3,size(Xq,1),size(Xq,2));
dM = dL; dN = dL;

% Try to do the same for rotated ndgrid points
[V, W] = differentiate(fitresult, Yrot, Zrot);
U = -ones(size(V));

% This part could potentially be optimized!
for ix=1:size(Xq,1)
   for iy=1:size(Xq,2)
      dN(:,ix,iy) = [U(ix,iy) V(ix,iy) W(ix,iy)];      
      dL(:,ix,iy) = cross(dN(:,ix,iy),z);
      dM(:,ix,iy) = cross(dL(:,ix,iy),dN(:,ix,iy));
      
      % Normalization
      dL(:,ix,iy) = dL(:,ix,iy) / norm(dL(:,ix,iy));
      dM(:,ix,iy) = dM(:,ix,iy) / norm(dM(:,ix,iy));
      dN(:,ix,iy) = dN(:,ix,iy) / norm(dN(:,ix,iy));
   end
end

%% Plots from rotated ndgrid in the local coordinate system

% box outputs
filename = '~/Ganymede/newPIC/run_G28_newPIC/box_FTE_G28_1200s.outs';

[~,~,fileinfo] = read_data(filename,'verbose',false);
npict = fileinfo.npictinfiles;
firstpict = 516; dpict = 10; lastpict = firstpict + 2*dpict; 

% Parameters and thresholds
Coef = 1.03; % expansion factor

% create new figure with specified size
hfig = figure('position',[10, 10, 800, 700]);
colormap(jet);

iplot = 1;
height = 0.24;

% Loop over snapshots
for ipict = firstpict:dpict:lastpict
   [filehead,data] = read_data(filename,'verbose',false,'npict',ipict);
   
   x = data.file1.x(:,:,:,1);
   y = data.file1.x(:,:,:,2);
   z = data.file1.x(:,:,:,3);
   
   bx = data.file1.w(:,:,:,5);
   by = data.file1.w(:,:,:,6);
   bz = data.file1.w(:,:,:,7);
   pe = data.file1.w(:,:,:,9);
   p  = data.file1.w(:,:,:,10);
   jx = data.file1.w(:,:,:,11);
   jy = data.file1.w(:,:,:,12);
   jz = data.file1.w(:,:,:,13);
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
   
   bxv= interp3(x, y, z, bx, Coef*Xq, Coef*Yq, Coef*Zq);
   byv= interp3(x, y, z, by, Coef*Xq, Coef*Yq, Coef*Zq); 
   bzv= interp3(x, y, z, bz, Coef*Xq, Coef*Yq, Coef*Zq);
   jxv= interp3(x, y, z, jx, Coef*Xq, Coef*Yq, Coef*Zq);
   jyv= interp3(x, y, z, jy, Coef*Xq, Coef*Yq, Coef*Zq); 
   jzv= interp3(x, y, z, jz, Coef*Xq, Coef*Yq, Coef*Zq);  
   pev= interp3(x, y, z, pe, Coef*Xq, Coef*Yq, Coef*Zq);
   pv = interp3(x, y, z, p,  Coef*Xq, Coef*Yq, Coef*Zq);
   
   jv = sqrt(jxv.^2 + jyv.^2 + jzv.^2);
      
   % Transform vectors into local coordinate system
   uL = Inf(size(Xq)); uM = uL; uN = uL;
   bL = Inf(size(Xq)); bM = bL; bN = bL;
      
   % This could potentially be improved!
   for ix=1:size(Xq,1)
      for iy=1:size(Xq,2)
         bL(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
         bM(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
         bN(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dN(:,ix,iy);
      end
   end
 
   % Bnormal plots
   h = subplot(4,3,iplot);
   
   contourf(Yq,Zq,bN,50,'Linestyle','none'); 
   title(sprintf('t=%ds',ipict));
   
   if iplot == 1
      text(1.3,0,'Bnormal [nT]','rotation',90,...
         'fontsize',14,'FontWeight','bold',...
         'horizontalalignment','center','verticalalignment','bottom');
   elseif iplot == 2
      h.YTickLabel = [];
   elseif iplot == 3
      h.YTickLabel = [];
      % Get the current axis size
      originalSize = get(gca, 'Position');
      colorbar; 
      set(gca,'Position', originalSize);
   end
   
   axis tight equal   
   h.Position(4) = height;
   h.XTickLabel = [];   
   h.XDir = 'reverse';
   h.FontSize = 12;
   caxis([-60 60])
   
   % Current density plots
   h = subplot(4,3,iplot+3);
   
   contourf(Yq,Zq,jv,50,'Linestyle','none'); 
   
   if iplot == 1
      text(1.3,0,'j [\mu A/m^2]','rotation',90,...
         'fontsize',14,'FontWeight','bold',...
         'horizontalalignment','center','verticalalignment','bottom');
   elseif iplot == 2
      h.YTickLabel = [];
      xlabel('Y [R_G]','FontWeight','bold','fontsize',14,...
         'horizontalalignment','center','verticalalignment','bottom');
      ylabel('Z [R_G]','FontWeight','bold','fontsize',14,...
         'horizontalalignment','center','verticalalignment','bottom');      
   elseif iplot == 3
      h.YTickLabel = [];
      % Get the current axis size
      originalSize = get(gca, 'Position');
      c = colorbar;
      set(gca,'Position', originalSize);
   end
   
   axis tight equal     
   h.Position(4) = height;
   h.XTickLabel = [];
   h.XDir = 'reverse';
   h.FontSize = 12;
   caxis([0.1 0.7])
   
   % Electron pressure plots
   h = subplot(4,3,iplot+6);
   
   contourf(Yq,Zq,pev,50,'Linestyle','none'); 

   if iplot == 1
      text(1.3,0,'P_e [nPa]','rotation',90,...
         'fontsize',14,'FontWeight','bold',...
         'horizontalalignment','center','verticalalignment','bottom');
   elseif iplot == 2
      h.YTickLabel = [];
   elseif iplot == 3
      h.YTickLabel = [];
      % Get the current axis size
      originalSize = get(gca, 'Position');
      colorbar;
      set(gca,'Position', originalSize);   
   end   

   axis tight equal   
   h.Position(4) = height;
   h.XDir = 'reverse';
   h.FontSize = 12;   
   caxis([0.2 4]);
     
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
h.YLim = [0.2,0.4];

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