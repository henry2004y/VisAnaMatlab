function galileo_b_plot( varargin )
%galileo_b_plot Processing Galileo observation/simulation data.
%   Read in Galileo trajectory and B data and plot the time series x and B.
%INPUT:
% flyby  : the index of flyby (default 8)
% Dir    : the directory of folder (default is current dir)
% DoMovie: logical for saving movies for trajectories (default 0)
% DoPlotB: logical for plotting magnetic observations (default 1)
% DoSaveB: logical for saving B plots to png          (default 0)
%
% Hongyang Zhou, hyzhou@umich.edu 06/13/2017 ver1.0
% Log
%----------------------------
% Modified on 1/1/2018 ver1.1
% 1. add Dir option to input arguments
%   
%--------------------------------------------------------------------------
   
% Input arguments
p = inputParser;
defaultDir     = '.';
defaultFlyby   = '8';
expectedFlybys = {'1','2','7','8','28','29'};
defaultDoMovie = 0;
defaultDoPlotB = 1;
defaultDoSaveB = 0;

addOptional(p,'flyby',defaultFlyby,...
   @(x) any(validatestring(x,expectedFlybys)));
addOptional(p,'DoMovie',defaultDoMovie,@isnumeric);
addOptional(p,'DoPlotB',defaultDoPlotB,@isnumeric);
addOptional(p,'DoSaveB',defaultDoSaveB,@isnumeric);
addOptional(p,'Dir',defaultDir,@ischar);

parse(p,varargin{:});


% Import position and mag. data
flybyfile=strcat('Galileo_G',p.Results.flyby,'_flyby_MAG.dat');
f = fullfile('~/Ganymede/GalileoData/galileomagdata',flybyfile);
delimiterIn = ' ';
headerlinesIn = 2;

G = importdata(f,delimiterIn,headerlinesIn);

Rg = 1;                     % Normalized Ganymede radius
N  = numel(G.data(:,1));    % Number of data
   
   %% Trajectory plot
   figure; hold on;  axis equal; 
   sphere; shading interp; material metal; 
   light('Style','infinite','Position',[-1 -2 2])
   h = plot3(G.data(:,7)/Rg,G.data(:,8)/Rg,G.data(:,9)/Rg);
   h.LineWidth = 1.2;
   hold off; grid on;
   xlabel('x'); ylabel('y'); zlabel('z');
   legend(h,{strcat('G',p.Results.flyby)},'Location','best');
   a = gca;
   a.LineWidth = 1.2; a.FontSize = 16;
   view(3)

   % Make avi movie
   if p.Results.DoMovie==1
      vobj=VideoWriter('orbit_Galileo', 'Motion JPEG AVI');
      vobj.FrameRate=40;
      vobj.Quality=75;
      open(vobj);
      for i=1:50:N
        figure(2);
        surf(x, y, z); axis equal;
        xlim([-5 8]); ylim([-10 10]);zlim([-5 5]);
        colormap winter; hold on
        for j=1:numfiles
           plot3(G{j}.data(:,7)./Rg,G{j}.data(:,8)./Rg,G{j}.data(:,9)./Rg,...
               'LineWidth',1.2);
        end
        xlabel('x'); ylabel('y'); zlabel('z');
        title('Galileo flyby');
        legend({'G1','G28','G29','G2','G7','G8'},'Location','best');
        for j=1:numfiles
           plot3(G{j}.data(i,7)./Rg,G{j}.data(i,8)./Rg,G{j}.data(i,9)./Rg,'o'); 
        end
        hold off
        F=getframe(2);
        writeVideo(vobj, F);
        cla(gca)
      end
      close(vobj)
   end

   %% B .vs. time
   if p.Results.DoPlotB==1
      % time in mins
      %time = dataFlyby(:,4)*60 + dataFlyby(:,5) + dataFlyby(:,6)/60;
      time = datetime(G.data(:,1:6));
      B    = G.data(:,10:12);
      figure('Position', [100, 100, 1049, 895]);
      subplot(411); LW1 = 3; LW2 = 1.5; FS=20;
      plot(time, B(:,1),'k','LineWidth',LW1); %Bx
      ylabel('Bx [nT]');
      title(strcat('Galileo G',p.Results.flyby,' Flyby Magnetic field'))
      set(gca,'FontSize',FS,'xMinorTick','on','LineWidth',LW2);
      subplot(412); 
      plot(time, B(:,2),'k','LineWidth',LW1); %By
      ylabel('By [nT]');
      set(gca,'FontSize',FS,'xMinorTick','on','LineWidth',LW2);
      subplot(413); 
      plot(time, B(:,3),'k','LineWidth',LW1); %Bz
      ylabel('Bz [nT]');
      set(gca,'FontSize',FS,'xMinorTick','on','LineWidth',LW2);
      subplot(414);  % |B|
      plot(time, sqrt(B(:,1).^2+B(:,2).^2+B(:,3).^2),'k','LineWidth',LW1); %|B|
      ylabel('B [nT]');
      set(gca,'FontSize',FS,'xMinorTick','on','LineWidth',LW2);

      if p.Results.DoSaveB==1
         OutputName = strcat('Galileo_G',p.Results.flyby,'_B.png');
         saveas(gcf,OutputName)
      end
   end
   
   
end

