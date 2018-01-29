function plot_data( varargin )
%PLOT_DATA Plot the variable from SWMF output.
%   INPUT
% data: original variable data
% filehead: header information
% vars: variable(s) for plotting
% plotmode:  (optional) type of plotting
% plotrange: (optional) range of plotting, [xmin xmax ymin ymax]
% plotinterval: (optional) interval for interpolation
% multifigure: logical value, 1 (default)for multifigure display,
%              0 for subplot
%
% Right now this can only deal with 2D plots. 3D streamline plots
% are demonstrated with fieldline.m and streamlinec.m
%
% 08/13
% I want to add a cut usage like in idl like
% cut = grid(10:30,*); plot_data
%
% 08/20
% I can add waitfor fcn for more control.
%
% 08/24
% Caution: there may be some issues with nx in generalized coords.!
%
% 08/29
% Think about how to simplify the string recognition for stream and quiver.
%
% 09/19
% Add an axis equal option.
%
% 09/22/2017
% For 3D visualization, I may make use of sliceomatic GUI. Try it!
%
% 01/19/2018
% For Cartesian Coordinates, I may use interp3 instead of
% scatteredInterpolant for speed. So the idea is to do interpolation onto
% all the points, and then choose the viewing region, instead of the
% inverse way.
% I need to switch the logics. Improve the code!
%
% 01/29/2018
% Right now this only works for one input data and header. Should I extend
% this to more?
%
% Hongyang Zhou, hyzhou@umich.edu ver1.0 06/14/2017
%--------------------------------------------------------------------------

%% Input parameters parser
p = inputParser;
defaultplotmode = 'contbar';
defaultplotrange = [];
defaultplotinterval = 1;
defaultmultifigure = 1;

addRequired(p,'data',@isstruct);
addRequired(p,'filehead',@isstruct);
addRequired(p,'func',@ischar);
addParameter(p,'plotrange',defaultplotrange,@isnumeric);
addParameter(p,'plotinterval',defaultplotinterval,@isnumeric);
addParameter(p,'plotmode',defaultplotmode,@ischar);
addParameter(p,'multifigure',defaultmultifigure,@islogical);

parse(p,varargin{:});

% I need to choose which variables to plot, as well as
% which plotmode to use. Think about it.
% func = strsplit(p.Results.func); % returns a cell array
% plotmode = strsplit(p.Results.plotmode); % returns a cell array
func = split(p.Results.func);
plotmode = split(p.Results.plotmode);
plotrange = p.Results.plotrange;
plotinterval = p.Results.plotinterval;
multifigure = p.Results.multifigure;
filehead = p.Results.filehead;
data = p.Results.data;
x = data.x; w = data.w;

%% Display parameters
disp('======= CURRENT PLOTTING PARAMETERS =======')

disp('======= PLOTTING PARAMETERS ===============')
fprintf('wnames =\n');
fprintf('%s ',filehead.wnames{:});
fprintf('\n======================================\n')

% Gabor uses this to do the plotting interactively
%read_plot_param

% Get transformation parameters and calculate grid
% Matlab cannot do contour plots in generalized coordinates, at least
% I don`t know how to do it. So I need a grid transformation whenever
% the data is in non-Cartesian coordinates. Just for now I use a if
% statement in the Plot loop to replace this.
%read_transform_param(filehead.gencoord)

% do_transform

%% Display min and max for each variable
funcs = strsplit(p.Results.func,{' ',';'});
for ifunc = 1:numel(funcs)
   VarIndex_ = strcmpi(funcs{ifunc},filehead.wnames);
   fprintf('Min and Max value for %s :%f, %f\n',funcs{ifunc}...
      ,min(w(:,1,VarIndex_)),max(w(:,1,VarIndex_)))
end

%% plot multiple variables with same plotmode 
if numel(plotmode) < numel(func)
   plotmode(end+1:numel(func)) = plotmode(1);
end

%% Plot
for ivar = 1:numel(func)
   % I need to think of a better way to check. now this cannot identify the
   % vars for streamline or quiver plotting!!!
   %validatestring(func(ivar),filehead.wnames)
   switch plotmode{ivar}
      case {'mesh','meshbar','meshbarlog','cont','contbar',...
            'contlog','contbarlog'}
         % find the index for var in filehead.wnames
         % I want to make this function more powerful to include
         % plotting derived variables, but it may not seem to be easy!
         VarIndex_ = strcmpi(func(ivar),filehead.wnames);
         W = reshape(w(:,:,VarIndex_),[],1);
         if isempty(W)
            error('%s not found in output variables!',func{ivar})
         end

         
         if multifigure; figure; end
         
         % Reorganize and pick data in plot region only
         % Somebody says it`s better to know the dimension ahead for the
         % sake of speed. Or maybe it`s even better not to do this
         % reshaping.
         X = reshape(x(:,:,1),[],1);
         Y = reshape(x(:,:,2),[],1);
         
         if ~isempty(plotrange)
            axis(plotrange)
         else
            plotrange(1) = min(X);
            plotrange(2) = max(X);
            plotrange(3) = min(Y);
            plotrange(4) = max(Y);
         end
         
         xyIndex = X > plotrange(1) & X < plotrange(2) & ...
            Y > plotrange(3) & Y < plotrange(4);
         X = X(xyIndex);
         Y = Y(xyIndex);
         W = W(xyIndex);
         
         if filehead.gencoord % Generalized coordinates
            % Default is linear interpolation
            F = scatteredInterpolant(X,Y,W);
            [xq,yq] = meshgrid(plotrange(1):plotinterval:plotrange(2)...
               ,plotrange(3):plotinterval:plotrange(4));
            vq = F(xq,yq);
            
            switch string(plotmode{ivar})
               case 'contbar'
                  contourf(xq,yq,vq,20,'Edgecolor','none'); c = colorbar;
                  %c.Label.String = '[nPa]';
                  %ylabel(c,'$\mu A/m^2$','Interpreter','LateX')
               case 'cont'
                  contour(xq,yq,vq,20,'Edgecolor','none');
               case 'contlog'
                  contourf(xq,yq,log(vq),20,'Edgecolor','none');
               case 'contbarlog'
                  contourf(xq,yq,log(vq),20,'Edgecolor','none');
                  c = colorbar;
                  c.Label.String = 'log10';
               case 'meshbar'
                  mesh(xq,yq,vq); colorbar
               case 'mesh'
                  mesh(xq,yq,vq);
               case 'meshbarlog'
                  mesh(xq,yq,log(vq)); c= colorbar;
                  c.Label.String = 'log10';
            end
         else % Cartesian coordinates, using meshgrid
            % Note: the number of points in each dimension can be
            % changed to the proper number you like.
            if isempty(plotrange)
               xlin = linspace(min(X),max(X),filehead.nx(1));
               ylin = linspace(min(Y),max(Y),filehead.nx(2));
            else
               xlin = linspace(plotrange(1),plotrange(2),filehead.nx(1));
               ylin = linspace(plotrange(3),plotrange(4),filehead.nx(2));
            end
            [xq,yq] = meshgrid(xlin,ylin);
            f = scatteredInterpolant(X,Y,W);
            vq = f(xq,yq);
            
            switch string(plotmode{ivar})
               case 'contbar'
                  contourf(xq,yq,vq,20,'Edgecolor','none'); c = colorbar;
               case 'cont'
                  contour(xq,yq,vq,'Edgecolor','none');
               case 'contlog'
                  contourf(xq,yq,log(vq),20,'Edgecolor','none');
               case 'contbarlog'
                  contourf(xq,yq,log(vq),20,'Edgecolor','none');
                  c = colorbar;
                  c.Label.String = 'log10';
               case 'meshbar'
                  mesh(xq,yq,vq); colorbar
               case 'mesh'
                  mesh(xq,yq,vq);
               case 'meshbarlog'
                  mesh(xq,yq,log(vq)); c= colorbar;
                  c.Label.String = 'log10';
            end
         end
         
         xlabel(filehead.variables{1}); ylabel(filehead.variables{2});
         title(filehead.wnames{VarIndex_});
         dim = [0.125 0.013 0.2 0.045];
         str = sprintf('it=%d, time=%.1fs',filehead.it,filehead.time);
         annotation('textbox',dim,'String',str,'FitBoxToText','on',...
            'FontWeight','bold');


      case {'trimesh','trisurf','tricont','tristream'} % triangular mesh
         figure;
         % find the index for var in filehead.wnames
         VarIndex_ = strcmpi(func(ivar),filehead.wnames);
         X = reshape(x(:,:,1),[],1);
         Y = reshape(x(:,:,2),[],1);
         W = reshape(w(:,:,VarIndex_),[],1);
         if ~isempty(plotrange)
            xyIndex = X > plotrange(1) & X < plotrange(2) ...
                    & Y > plotrange(3) & Y < plotrange(4);           
            X = X(xyIndex);
            Y = Y(xyIndex);
            W = W(xyIndex);
            
%             %tristream test
%             u = reshape(w(:,:,2),1,[]);
%             v = reshape(w(:,:,4),1,[]);
%             u = u(xyIndex);
%             v = v(xyIndex);
         end
         t = delaunayn([X, Y]);
         if string(plotmode{ivar}) == 'trimesh'
            trimesh(t,X,Y,W);
         elseif string(plotmode{ivar}) == 'trisurf'
            trisurf(t,X,Y,W,'EdgeColor','none');
         elseif string(plotmode{ivar}) == 'tristream'
%             x0 = ones(7,1)*-3;
%             y0 = [ -3 -2 -1 0 1 2 3]';
%             T = triangulation(t,x,y);
%             %triplot(T);
%             %FlowP=TriStream(T,u,v,x0,y0,1,2e3);
%             FlowP=TriStream(T,u,v,x0,y0);
%             PlotTriStream(FlowP,'r');
         else
            % I need to use patch to write my own tricont function!
            %tricont(t,x,y,w);
         end
         
      case {'stream','streamover'}
         % stream is not working!!! Fix it!
         if strcmp(plotmode{ivar},'streamover')
            % Overplotting with more variables, keyword 'over'
            hold on
         else
            if multifigure; figure; end
         end
         
         % find the index for var in filehead.wnames
         VarStream  = split(func(ivar),';');
         VarIndexS1 = strcmpi(VarStream(1),filehead.wnames);
         VarIndexS2 = strcmpi(VarStream(2),filehead.wnames);
         
         F1 = scatteredInterpolant(x(:,1,1),x(:,1,2),w(:,1,VarIndexS1));
         v1 = F1(xq,yq);
         F2 = scatteredInterpolant(x(:,1,1),x(:,1,2),w(:,1,VarIndexS2));
         v2 = F2(xq,yq);
         
         s = streamslice(xq,yq,v1,v2);
         
         for is=1:numel(s)
            s(is).Color = 'w'; % Change streamline color to white
            s(is).LineWidth = 1.5;
         end
                  
         if ~isempty(plotrange) ;axis(plotrange); end
         
      case {'quiver','quiverover'}
         if strcmp(plotmode{ivar},'quiverover')
            % Overplotting with more variables, keyword 'over'
            hold on
         else
            if multifigure; figure; end
         end         
         
         % find the index for var in filehead.wnames
         VarQuiver  = split(func(ivar),';');
         VarIndexS1 = strcmpi(VarQuiver(1),filehead.wnames);
         VarIndexS2 = strcmpi(VarQuiver(2),filehead.wnames);
         
         q = quiver(x(:,1,1),x(:,1,2),w(:,1,VarIndexS1),w(:,1,VarIndexS2));
         %q.Color = 'w';
         q.AutoScaleFactor = 0.5;
         if ~isempty(plotrange); axis(plotrange); end
         
      case 'grid'    % Grid plot
         figure
         scatter(x(:,1,1),x(:,1,2),'.','LineWidth',0.1);
         
         if ~isempty(plotrange) 
            axis(plotrange)
         end

         xlabel(filehead.variables{1}); ylabel(filehead.variables{2});
         title('Grid illustration');
         dim = [0.125 0.013 0.1 0.046];
         str = sprintf('it=%d',filehead.it);
         annotation('textbox',dim,'String',str,'FitBoxToText','on');
         %drawnow     
      otherwise
         error('unknown plot mode: %s',plotmode{ivar})
   end
 
   % There`s a problem with this two lines: if there are several figures, 
   % then only the last one has the expected font size! Originally this was
   % at the bottom for loop end, but I moved it upward inside the inner
   % loop.
   ax = gca;
   ax.FontSize = 16;  
end


end
