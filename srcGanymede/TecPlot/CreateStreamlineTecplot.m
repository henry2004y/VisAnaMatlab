function CreateStreamlineTecplot(varargin)
% Matlab function for generating macro to be used in Tecplot for creating
% B-field lines.
%
% Note that the user needs to change Var number according to their outputs.
% The defaults are
% Bx: Var number 8
% By: Var number 9
% Bz: Var number 10
% Change these according to your dataset.
%
%INPUT ARGUMENT:
% filename: filename for .mcr
% filetype: file types for seed picking (2D/3D/line/boundarycurve)
% radius: radial location to start streamlines in R_planet
% bidirection: logical for field line tracing direction.
%
%OUTPUT
% macro in Tecplot format
% 
%
% Example 1: creating field lines at fixed latitudes/longitudes over a
% sphere in 3D
% CreateStreamlinesTecplot('filename','CreateStreamlines.mcr',...
%    'filetype','3D','radius',2)
%
% Example 2: creating field lines across a line of point
% CreateStreamlinesTecplot('filetype','line','xStart',[-3 0 0],...
%    'xEnd',[-2 0 0])
%
% Example 3: creating field lines in 2D (for y=0 XZ plane only)
%
%
%
% Happy plotting! :)
% Written by Yash Sarkango,
% Modified by Hongyang Zhou, hyzhou@umich.edu 02/12/2018
%--------------------------------------------------------------------------

%  if nargin==0
%     [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
%     v = x .* exp(-x.^2 - y.^2 - z.^2);
%     sliceomatic(v)
%     return
%  end
 
%% Input parameters parser
p = inputParser;
defaultfilename = 'CreateStreamlines.mcr';
defaultfiletype = '3D'; 
expectedtypes = {'3D','3d','2D','2d','line','boundarycurve'};
defaultradius = 1.0;
defaultbidirection = true;

%addRequired(p,'filename',defaultfilename,@ischar);
addParameter(p,'filename',defaultfilename,@ischar);
%addOptional(p,'npict',defaultnpict,@isnumeric);

addParameter(p,'filetype',defaultfiletype,...
   @(x) any(validatestring(x,expectedtypes)));
addOptional(p,'radius',defaultradius,@isnumeric);
addOptional(p,'xStart',[0 0 0],@(x)numel(x)==3);
addOptional(p,'xEnd',[0 0 0],@(x)numel(x)==3);
addOptional(p,'nInterval',10,@isnumeric);
addOptional(p,'curvetype','day',@ischar);
addParameter(p,'direction',defaultbidirection,@islogical);

parse(p,varargin{:});
p = p.Results;

%%

% Specify file name and open for writing
File.name    = p.filename;
File.ID      = fopen(File.name,'w');

 
switch p.filetype
    % ----- 3D Streamlines ------------------------------------------------  
    case {'3D','3d'}
        N_phi               = 12;           % No. of points in azimuth
        phi                 = 0:2*pi/N_phi:(2*pi - 2*pi/N_phi);

        %theta               = [20,30,40,50,60,65,70,75,77,78,79.9,80];
        %theta = [20,40,60,80,140,150];
        theta = [30,45,70,115,123,160];
        
        fprintf(File.ID,'#!MC 1410\n');
        fprintf(File.ID,'$!VarSet |MFBD| = ''/''\n');
        fprintf(File.ID,'$!GLOBALTHREEDVECTOR UVAR = 8\n');
        fprintf(File.ID,'$!GLOBALTHREEDVECTOR VVAR = 9\n');
        fprintf(File.ID,'$!GLOBALTHREEDVECTOR WVAR = 10\n');
        fprintf(File.ID,'$!STREAMTRACELAYERS SHOW = YES\n');
        fprintf(File.ID,'$!FRAMESETUP DEFAULT3DSTREAMTRACE{STREAMTYPE = VOLUMELINE}\n');
        fprintf(File.ID,'$!STREAMTRACE RESETDELTATIME\n');
        fprintf(File.ID,'$!STREAMATTRIBUTES RODRIBBON{WIDTH = 32.9090909090908923}\n');

        for i = 1:max(size(theta))
           for j = 1:max(size(phi))
              x    = p.radius*sind(theta(i))*cos(phi(j));
              y    = p.radius*sind(theta(i))*sin(phi(j));
              z    = p.radius*cosd(theta(i));
              fprintf(File.ID,'$!STREAMTRACE ADD\n');
              fprintf(File.ID,'  STREAMTYPE = VOLUMELINE\n');
              fprintf(File.ID,'  STREAMDIRECTION = BOTH\n');
              fprintf(File.ID,'  STARTPOS\n');
              fprintf(File.ID,'  {\n');
              fprintf(File.ID,'  X = %.10f\n',x);
              fprintf(File.ID,'  Y = %.10f\n',y);
              fprintf(File.ID,'  Z = %.10f\n',z);
              fprintf(File.ID,'  }\n');
           end
        end

    % ----- 2D Streamlines ------------------------------------------------    
    case {'2D','2d'}
        
        phi                 = [0 pi];       % 2 phi values in XZ plane

        theta               = pi/180*linspace(0,90,30); % change as needed
        
        fprintf(File.ID,'#!MC 1410\n');
        fprintf(File.ID,'$!VarSet |MFBD| = ''/''\n');
        fprintf(File.ID,'$!GLOBALTWODVECTOR UVAR = 8\n');
        fprintf(File.ID,'$!GLOBALTWODVECTOR VVAR = 10\n');
        fprintf(File.ID,'$!STREAMTRACELAYERS SHOW = YES\n');
        fprintf(File.ID,'$!STREAMTRACE RESETDELTATIME\n');

        for i = 1:max(size(theta))
           for j = 1:max(size(phi))
              x    = r*sin(theta(i))*cos(phi(j));
              y    = r*sin(theta(i))*sin(phi(j));
              z    = r*cos(theta(i));
              fprintf(File.ID,'$!STREAMTRACE ADD\n');
              fprintf(File.ID,'  STREAMTYPE = TWODLINE\n');
              fprintf(File.ID,'  STREAMDIRECTION = BOTH\n');
              fprintf(File.ID,'  STARTPOS\n');
              fprintf(File.ID,'  {\n');
              fprintf(File.ID,'  X = %.10f\n',x);
              fprintf(File.ID,'  Y = %.10f\n',z);     %Y is 'z' in XZ plane
              fprintf(File.ID,'  }\n');
           end
        end
        
   % 3D seeds on curves ---------------------------------------------------
   case {'boundarycurve'}
      
      filename= '../newPIC/run_G8_newPIC/box_CPCP_G8_1200s.outs';
      ipict = 1;
      xCenter = 0; yCenter = 0; % coordinates of the center of points
      fprintf('ipict=%d\n',ipict);
      [~,data] = read_data(filename,'verbose',false,'npict',ipict);
      
      data = data.file1;
      x = data.x(:,:,:,1);
      y = data.x(:,:,:,2);
      status = data.w(:,:,:,14);
      status(status>1) = 2;
      status(status<1) = 0;
      
      x = x(status==2);
      y = y(status==2);
      
      % Find the boundary of magnetosphere
      k = boundary(x,y,1);
      % Integrate along the boundary
      x = x(k);
      y = y(k);
    
      if curvetype == 'day'
         % Pick the upstream half curve (use the line x==1 as boundary)
         k = x <= 1;       
         x = x(k);
         y = y(k);       
         
         % Re-order of boundary points in a counter-clockwise way
         angles = atan2d( y-yCenter,x-xCenter );
         sortIndexes = angles<0;
         angles(sortIndexes) = angles(sortIndexes) + 360;
         [sortAngles, sortIndexes] = sort(angles);
         x = x(sortIndexes);
         y = y(sortIndexes);
         z(1:numel(x)) = 2;
      elseif curvetype == 'night'
         % Pick the downstream half curve
         k = x > 1.1; % just for a test
         x = x(k);
         y = y(k);        
      end
      
      fprintf(File.ID,'#!MC 1410\n');
      fprintf(File.ID,'$!VarSet |MFBD| = ''/''\n');
      fprintf(File.ID,'$!GLOBALTHREEDVECTOR UVAR = 8\n');
      fprintf(File.ID,'$!GLOBALTHREEDVECTOR VVAR = 9\n');
      fprintf(File.ID,'$!GLOBALTHREEDVECTOR WVAR = 10\n');
      fprintf(File.ID,'$!STREAMTRACELAYERS SHOW = YES\n');
      fprintf(File.ID,'$!FRAMESETUP DEFAULT3DSTREAMTRACE{STREAMTYPE = VOLUMELINE}\n');
      fprintf(File.ID,'$!STREAMTRACE RESETDELTATIME\n');
      fprintf(File.ID,'$!STREAMATTRIBUTES RODRIBBON{WIDTH = 32.9090909090908923}\n');
      
      for i = 1:10:numel(x)
         fprintf(File.ID,'$!STREAMTRACE ADD\n');
         fprintf(File.ID,'  STREAMTYPE = VOLUMELINE\n');
         fprintf(File.ID,'  STREAMDIRECTION = BOTH\n');
         fprintf(File.ID,'  STARTPOS\n');
         fprintf(File.ID,'  {\n');
         fprintf(File.ID,'  X = %.10f\n',x(i));
         fprintf(File.ID,'  Y = %.10f\n',y(i));
         fprintf(File.ID,'  Z = %.10f\n',z(i));
         fprintf(File.ID,'  }\n');
      end
      
   % seeds along a line ---------------------------------------------------
   case 'line'
      
      x=linspace(p.xStart(1),p.xEnd(1),p.nInterval);
      y=linspace(p.xStart(2),p.xEnd(2),p.nInterval);
      z=linspace(p.xStart(3),p.xEnd(3),p.nInterval);

      fprintf(File.ID,'#!MC 1410\n');
      fprintf(File.ID,'$!VarSet |MFBD| = ''/''\n');
      fprintf(File.ID,'$!GLOBALTHREEDVECTOR UVAR = 8\n');
      fprintf(File.ID,'$!GLOBALTHREEDVECTOR VVAR = 9\n');
      fprintf(File.ID,'$!GLOBALTHREEDVECTOR WVAR = 10\n');
      fprintf(File.ID,'$!STREAMTRACELAYERS SHOW = YES\n');
      fprintf(File.ID,'$!FRAMESETUP DEFAULT3DSTREAMTRACE{STREAMTYPE = VOLUMELINE}\n');
      fprintf(File.ID,'$!STREAMTRACE RESETDELTATIME\n');
      fprintf(File.ID,'$!STREAMATTRIBUTES RODRIBBON{WIDTH = 32.9090909090908923}\n');

      for i = 1:p.nInterval
         fprintf(File.ID,'$!STREAMTRACE ADD\n');
         fprintf(File.ID,'  STREAMTYPE = VOLUMELINE\n');
         fprintf(File.ID,'  STREAMDIRECTION = BOTH\n');
         fprintf(File.ID,'  STARTPOS\n');
         fprintf(File.ID,'  {\n');
         fprintf(File.ID,'  X = %.10f\n',x(i));
         fprintf(File.ID,'  Y = %.10f\n',y(i));
         fprintf(File.ID,'  Z = %.10f\n',z(i));
         fprintf(File.ID,'  }\n');
      end

   otherwise
      error('Wrong filetype!')
end

fprintf('...Done.\n');

end