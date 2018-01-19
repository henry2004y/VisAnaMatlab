function plot_log_data( varargin )
%plot_log_data Plot the information from log file
%   INPUT
% data: original variable data
% filehead: header information
% vars: variable(s) for plotting
% plotmode:  (optional) type of plotting {'line','scatter'}
% plotrange: (optional) range of plotting
%
%--------------------------------------------------------------------------

% Input parameters parser
p = inputParser;
defaultplotmode = 'line';
defaultplotrange = [];

addRequired(p,'data',@isnumeric);
addRequired(p,'filehead',@isstruct);
addRequired(p,'func',@ischar);
addParameter(p,'plotrange',defaultplotrange,@isnumeric);
addParameter(p,'plotmode',defaultplotmode,@ischar);

parse(p,varargin{:});

func      = strsplit(p.Results.func);
plotmode  = strsplit(p.Results.plotmode);
plotrange = p.Results.plotrange;
filehead  = p.Results.filehead;
data      = p.Results.data;


for ivar = 1:numel(func)
   switch string(plotmode{ivar})
      case 'line'
         % find the index for var in filehead.variables
         VarIndex = strcmpi(func{ivar},filehead.variables);
         
         if all(VarIndex==0)
            error('unknown plotting variable.')
         end
         
         figure;
         
         plot(data(:,1),data(:,VarIndex))
         
         xlabel(filehead.variables{1});
         ylabel(filehead.variables{VarIndex});
         title('log file data');
         set(gca,'FontSize',14);
      case 'scatter'
         % find the index for var in filehead.variables
         VarIndex = strcmpi(func{ivar},filehead.variables);

         if all(VarIndex==0)
            error('unknown plotting variable.')
         end
         
         figure;
         
         scatter(data(:,1),data(:,VarIndex))
         
         xlabel(filehead.variables{1});
         ylabel(filehead.variables{VarIndex});
         title('log file data');
         set(gca,'FontSize',14);
         
      otherwise
         error('unknown plot mode for plot_log_data!')
   end
      
end

end