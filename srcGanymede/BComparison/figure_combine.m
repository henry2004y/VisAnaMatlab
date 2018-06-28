% Script for SWMF Galileo Flyby Simulation Comparison
% Purpose: Plot quasi-steady state magnetic field comparisons for all six
%          Galileo flybys.
%
% Prerequisite: Folder FileIO_SWMF must be included into Matlab path
%
% Hope to get a nice combines figure.
%
% Hongyang Zhou, hyzhou@umich.edu 02/06/2018

clear;clc
%% Parameters
flyby = [1 2 7 8 28 29];
height = 0.053; % height for each subplot

%%
figure('Unit','inches','Position', [0, 0, 7, 10]);
%figure('PaperType','usletter')
%figure('Position',[0,0,1400,2000])
LW1 = 1.5; LW2 = 1.5; FS=8;

% G1
h(1) = subplot(14,2,1);
h(2) = subplot(14,2,3);
h(3) = subplot(14,2,5);
h(4) = subplot(14,2,7);


% G2
h(5) = subplot(14,2,2);
h(6) = subplot(14,2,4);
h(7) = subplot(14,2,6);
h(8) = subplot(14,2,8);

% G7
h(9) = subplot(14,2,11);
h(10) = subplot(14,2,13);
h(11) = subplot(14,2,15);
h(12) = subplot(14,2,17);

% G8
h(13) = subplot(14,2,12);
h(14) = subplot(14,2,14);
h(15) = subplot(14,2,16);
h(16) = subplot(14,2,18);

% G28
h(17) = subplot(14,2,21);
h(18) = subplot(14,2,23);
h(19) = subplot(14,2,25);
h(20) = subplot(14,2,27);

% G29
h(21) = subplot(14,2,22);
h(22) = subplot(14,2,24);
h(23) = subplot(14,2,26);
h(24) = subplot(14,2,28);


%%
for iflyby = 1:6
   % Read observation data
   flybyfile = strcat('Galileo_G',int2str(flyby(iflyby)),'_flyby_MAG.dat');
   f = fullfile('~/Ganymede/GalileoData/galileomagdata',flybyfile);
   [~,data] = read_log_data(f);
   
   time = datetime(data(:,1:6));
   xyz  = data(:,7:9);
   Bobs = data(:,10:12);
   BobsStrength = sqrt(Bobs(:,1).^2+Bobs(:,2).^2+Bobs(:,3).^2);
   
   % Read simulation data
   filename = strcat('~/Ganymede/MOP2018/runG',int2str(flyby(iflyby)),...
      '_steady/GM/box*');
   [~,~,info] = read_data(filename);
   [filehead,data] = read_data(filename,'npict',info.npictinfiles,'verbose',false);
   
   % Interpolate simulation data to observation data
   data = data.file1;
   nx = filehead.nx; nw = filehead.nw;
   
   x = data.x(:,:,:,1);
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);
   
   w = data.w;
   % Find the correct index of variables
   bx_ = strcmpi('bx',filehead.wnames);
   by_ = strcmpi('by',filehead.wnames);
   bz_ = strcmpi('bz',filehead.wnames);
   bx = data.w(:,:,:,bx_);
   by = data.w(:,:,:,by_);
   bz = data.w(:,:,:,bz_);
   
   % From ndgrid to meshgrid format
   x  = permute(x,[2 1 3]);
   y  = permute(y,[2 1 3]);
   z  = permute(z,[2 1 3]);
   bx = permute(bx,[2 1 3]);
   by = permute(by,[2 1 3]);
   bz = permute(bz,[2 1 3]);
   
   % gridded interpolation, much faster than scattered interpolation
   Bsim = Inf(size(xyz,1),3);
   Bsim(:,1) = interp3(x,y,z,bx,xyz(:,1),xyz(:,2),xyz(:,3));
   Bsim(:,2) = interp3(x,y,z,by,xyz(:,1),xyz(:,2),xyz(:,3));
   Bsim(:,3) = interp3(x,y,z,bz,xyz(:,1),xyz(:,2),xyz(:,3));
   
   
   BsimStrength = sqrt(Bsim(:,1).^2+Bsim(:,2).^2+Bsim(:,3).^2);

   
   % Plot
   plot(h(1+(iflyby-1)*4),time, Bobs(:,1),'k',time,Bsim(:,1),...
      'r','LineWidth',LW1); %Bx
   plot(h(2+(iflyby-1)*4),time, Bobs(:,2),'k',time,Bsim(:,2),...
      'r','LineWidth',LW1); %By
   plot(h(3+(iflyby-1)*4),time, Bobs(:,3),'k',time, Bsim(:,3),...
      'r','LineWidth',LW1); %Bz
   plot(h(4+(iflyby-1)*4),time,BobsStrength,'k',time,BsimStrength,...
      'r','LineWidth',LW1); %|B|
   
   h(1+(iflyby-1)*4).XLim = [min(time) max(time)];
   h(2+(iflyby-1)*4).XLim = [min(time) max(time)];
   h(3+(iflyby-1)*4).XLim = [min(time) max(time)];
   h(4+(iflyby-1)*4).XLim = [min(time) max(time)]; 
   h(1+(iflyby-1)*4).YLim = [min(Bobs(:,1))-50 max(Bobs(:,1))+50 ];
   h(2+(iflyby-1)*4).YLim = [min(Bobs(:,2))-50 max(Bobs(:,2))+50 ];
   h(3+(iflyby-1)*4).YLim = [min(Bobs(:,3))-50 max(Bobs(:,3))+50 ];
   h(4+(iflyby-1)*4).YLim = [min(BobsStrength)-50 max(BobsStrength)+50 ];
end

for ih=1:24
   h(ih).FontSize = FS;
   h(ih).XAxis.MinorTick = 'on';
   h(ih).YMinorTick = 'on';
   h(ih).LineWidth = LW2;
   
   switch ih
      case {1,5,9,13,17,21}
         h(ih).YLabel.String = 'Bx [nT]';
         h(ih).XTickLabel = [];
         h(ih).Position(4) = height;
         %h(ih).XLim = [min(time) max(time)];
         h(ih).Title.String = {strcat('(',char(int8(ih/4)+97),')');...
            strcat('Galileo G',int2str(flyby(int8(ih/4)+1)),...
            ' Flyby Magnetic field')};
         if ih==1
            legend(h(ih),{'Obs','Sim'})
            h(ih).Legend.Position = [0.457 0.95 0.108 0.031];
         end
      case {2,6,10,14,18,22}
         h(ih).YLabel.String = 'By [nT]';
         h(ih).XTickLabel = [];
         h(ih).Position(4) = height;
      case {3,7,11,15,19,23}
         h(ih).YLabel.String = 'Bz [nT]';
         h(ih).XTickLabel = [];
         h(ih).Position(4) = height;
      case {4,8,12,16,20,24}
         h(ih).YLabel.String = 'B [nT]';
         h(ih).Position(4) = height;
   end
end

%print('-bestfit','BestFitFigure','-dpdf')
saveas(gcf,'BField','epsc')
%print('test','-depsc')
%print('test2','-dtiff')
%print('test3','-dpdf')