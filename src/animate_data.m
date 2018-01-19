function animate_data( varargin )
%animate_data Plot picture(s) as animation for SWMF output files.  
%   INPUT
% filename: file name(s)
% func: variable(s) for animation
% plotmode:  (optional) type of plotting
% plotrange: (optional) range of plotting, [xmin xmax ymin ymax]
% plotinterval: (optional) interval for interpolation
% firstpict: first picture to be plotted
% dpict: interval of pictures to be plotted
% npictmax: max number of pictures to be plotted
% savemovie: logical for saving, 1 save, 0 no save
% movieformat: string of saving format, default is 'avi'
%
% I am planning to use javascript for visualization, and my own script
% for saving into movies. But I need to figure out some better way to make
% movies with good quality and small size.
%
% Another problem is the uniform colorbar for animation. I need to fix
% this.
%
% set(gca, 'nextplot', 'replacechildren');
% caxis manual;  % allow subsequent plots to use the same color limits
% caxis([-1 1]); % set the color axis scaling to your min and max color limits
%
% Right now this can only deal with 2D plots, only one output file!
% Another problem is the colorbar: I should use a uniform colorbar for
% every snapshot in the animation!
%
% Hongyang Zhou, hyzhou@umich.edu
%--------------------------------------------------------------------------

% Input parameters parser
p = inputParser;
defaultplotmode = 'contbar';
defaultplotrange = [];
defaultplotinterval = 1;
defaultdpict = 1;
defaultfirstpict = 1;
defaultnpictmax = 100;
defaultsavemovie = false;
defaultmovieformat = 'avi';
expectedformat = {'avi','mp4','mj2'};

addRequired(p,'filename',@ischar);
addRequired(p,'func',@ischar);
addParameter(p,'plotrange',defaultplotrange,@isnumeric);
addParameter(p,'plotinterval',defaultplotinterval,@isnumeric);
addParameter(p,'plotmode',defaultplotmode,@ischar);
addParameter(p,'dpict',defaultdpict,@isnumeric);
addParameter(p,'firstpict',defaultfirstpict,@isnumeric);
addParameter(p,'npictmax',defaultnpictmax,@isnumeric);
addParameter(p,'savemovie',defaultsavemovie,@islogical);
addParameter(p,'movieformat',defaultmovieformat,...
   @(x) any(validatestring(x,expectedformat)));

parse(p,varargin{:});

filename    = p.Results.filename;
func        = p.Results.func;
plotmode    = p.Results.plotmode;
plotrange   = p.Results.plotrange;
plotinterval= p.Results.plotinterval;
dpict       = p.Results.dpict;
firstpict   = p.Results.firstpict;
npictmax    = p.Results.npictmax;
savemovie   = p.Results.savemovie;
movieformat = p.Results.movieformat;

% Get info of file
[~,~,filelist] = read_data(filename);

npict = floor( (filelist.npictinfiles-firstpict)/dpict+1 );
if npict > npictmax; npict = npictmax; end
if npict < 0; npict = 0; end

if npict==0 
   error(['There are no frames to animate! '...
      'Check the following settings:\n'...
      '   npictinfiles = %d\n'...
      '   firstpict    = %d\n'...
      '   dpict        = %d\n'...
      '   npictmax     = %d\n'],...
      filelist.npictinfiles,firstpict,dpict,npictmax)
end

ipict = 1; icount = 1;
if savemovie
   !mkdir -p Movie
   switch string(movieformat)
      case 'avi'
         v = VideoWriter('Movie/movie.avi');
      case 'mp4'
         v = VideoWriter('Movie/movie.mp4');
      case 'mj2'
         v = VideoWriter('Movie/movie.mj2');         
   end   
   v.FrameRate = 10;
   v.open
   set(gca,'nextplot','replacechildren');
   
   while icount <= npict
      [filehead,data] = read_data(filename,'npict',ipict,'verbose',false);
      plot_data(data.file1,filehead,func...
         ,'plotmode',plotmode,'plotrange',plotrange...
         ,'plotinterval',plotinterval,'multifigure',false);
      frame = getframe(gcf);
      
      %clf;
      writeVideo(v,frame);     
      ipict = ipict + dpict; icount = icount + 1;
   end   
   v.close
   close all
else
   F(npict) = struct('cdata',[],'colormap',[]);
   while icount <= npict
      [filehead,data] = read_data(filename,'npict',ipict,'verbose',false);
      plot_data(data.file1,filehead,func...
         ,'plotmode',plotmode,'plotrange',plotrange...
         ,'plotinterval',plotinterval,'multifigure',false);
      F(icount) = getframe(gcf);
      %set(gca,'nextplot','replacechildren');
      clf;      
      ipict = ipict + dpict; icount = icount + 1;
   end
   % I don`t know why, but the display is shifted?
   movie(F,1)
end

end
