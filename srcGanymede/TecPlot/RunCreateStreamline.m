% Run script for creating streamlines for Ganymede G8 3d outputs
%
%

% Upstream half curve
CreateStreamlineTecplot('filename','CreateStreamlines_day.mcr',...
   'filetype','boundarycurve')

% Downstream half curve
CreateStreamlineTecplot('filename','CreateStreamlines_night.mcr',...
   'filetype','boundarycurve','curvetype','night')

% lines
CreateStreamlineTecplot('filename','CreateStreamlines_line.mcr',...
   'filetype','line','xStart',[-4 0 0],'xEnd',[-2.5 0 0],'nInterval',10)

%%
CreateStreamlineTecplot('filename','TecPlot/CreateStreamlines_line.mcr',...
   'filetype','line','xStart',[-2.8 0 0],'xEnd',[-1.5 0 0],'nInterval',13)

CreateStreamlineTecplot('filename','TecPlot/CreateStreamlines_line2.mcr',...
   'filetype','line','xStart',[2 0 0],'xEnd',[5 0 0],'nInterval',10)


%
CreateStreamlineTecplot('filename','TecPlot/CreateStreamlines_sphere.mcr',...
   'filetype','3d')

%% 
CreateStreamlineTecplot('filename','TecPlot/CreateStreamlines_sphere_G28.mcr',...
   'filetype','3d')