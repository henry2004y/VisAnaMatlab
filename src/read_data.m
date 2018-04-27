function [ filehead,data,filelist ] = read_data( filename, varargin )
%read_data Read data from BATSRUS output files
%   Read the npict-th snapshot from an ascii or binary data file into          
%   the x (coordinates) and w (data) arrays.
%
% Usage:                                                                     
%                                                                               
% filename='...' ; set file to read from (optional)                             
% npict=...      ; set snapshot index (optional), default 1                            
% read_data
%
%   The "x" and "w" arrays and the header info will be read from the file.
%
%   The same npict-th snapshot can be read from up to 10 files by e.g.
%   setting                                                                       
% filename='data/file1.ini data/file2.out'
%
%   In this case the data is read into x0,w0 and x1,w1 for the two files
%--------------------------------------------------------------------------

% Right now I cannot deal with multiple files. 06/25
% Do it later!!!
% I am confused with nx: what is it? In my current application it is always
% 1, but I may need to change it in the future for more general files.
% 07/03
%
% The filetype input is always overwritten inside the function,
% which means that it is useless to give it as an argument!
% Consider remove it!
% 07/31
%
% Removed the global variables. 08/19/2017
%
% Removed preallocation of filehead in the main function: move to a
% get_file_head as its persistent variable. 08/20/2017  
%
% I decided to use a slightly different reading logic compared to Gabor`s.
% 08/21/2017
%
% A more audacious trial: use nested function to save memory. I will create
% another branch for this, just for testing. 08/27/2017
%
% Now I understand why shouldn`t I use the first header! iteration number
% and time are different! 08/29/2017
%
% I can try containers.map or table data structure to save vars. 09/22/2017
%
% I need a way to deal with no output arguments. 02/14/2018
%
%--------------------------------------------------------------------------


%% Input parameters parser
p = inputParser;
%defaultfilenames = '';       % array of filenames
%defaultnfile = 0;            % number of files
defaultfiletype = 'binary';  % file types
defaultnpict = 1;            % index of snapshot to be read
defaultnpictinfiles = 1;     % number of pictures in each file
expectedtypes = {'binary','real4','ascii','log'};
defaultverbose = true;

addRequired(p,'filenames',@ischar);
addOptional(p,'npict',defaultnpict,@isnumeric);

addParameter(p,'filetype',defaultfiletype,...
   @(x) any(validatestring(x,expectedtypes)));
addParameter(p,'npictinfiles',defaultnpictinfiles,@isnumeric);
addParameter(p,'verbose',defaultverbose,@islogical);

parse(p,filename,varargin{:});


%% Check the existence of files
%
% In this way I cannot have several filenames separated with space
% First do one file, then try several files.
% There`s a bug for using * here!
filename = split(filename);
%filelist = struct('name',{});
for ifile=1:size(filename,1)
   newfile = dir(char(filename(ifile,:)));
   if isempty(newfile)
      error('Error in read_data: no matching filename was found for %s.'...
         ,filename{ifile})
   end
   for inewfile=1:numel(newfile)
      try
         filelist(end+1) = newfile(inewfile);
      catch
         filelist = newfile(inewfile);
      end
   end
end

nfile = numel(filelist);
%filehead(nfile).headline = '';
if nfile>2
   fprintf('nfile = %d\n',nfile);
   fprintf('filenames = %s\n',filelist.name);
   error('Error in read_data: cannot handle more than 2 files.')
end

[filelist,fileID,pictsize] = ...
   get_file_types(nfile,filelist);

if p.Results.verbose
   fprintf('filename=%s\n',filelist.name);
   fprintf('npict=%d\n',filelist.npictinfiles);
end

for ifile=1:nfile
   if any(bsxfun(@minus, filelist(ifile).npictinfiles, p.Results.npict)<0)
      error('file %d: npict out of range!',ifile)
   end
   frewind(fileID(ifile));
end

%% Read data from file ifile
for ifile=1:nfile
   if strcmp(filelist(ifile).type,'log')     
      [filehead,data] = read_log_data(fullfile(filelist(ifile).folder,...
         filelist(ifile).name));     
      return
   else          
      % Skip npict-1 snapshots (because we only want npict snapshot)      
      fseek(fileID(ifile),...
         pictsize(ifile)*(p.Results.npict-1),'cof');  

      [~,filehead(ifile)] = get_file_head(fileID(ifile),...
         filelist(ifile).type);
      
      % Read data
      switch string(lower(filelist(ifile).type))
         case 'ascii'
            [x,w] = get_pict_asc ( fileID(ifile),filehead(ifile) );
         case 'binary'
            [x,w] = get_pict_bin ( fileID(ifile),filehead(ifile) );
         case 'real4'
            [x,w] = get_pict_real( fileID(ifile),filehead(ifile) );
         case 'log'
         otherwise
            error('get_pict: unknown filetype: %s',filelist(ifile).type)
      end
      
      set_units(filehead(ifile),' ');
      
      if p.Results.verbose
         show_head(filelist(ifile),ifile,filehead(ifile))
      end
      
      % Produce a wnames from the last file
      filehead(ifile).wnames=filehead(ifile).variables(...
         filehead(ifile).ndim+1:filehead(ifile).ndim+filehead(ifile).nw);

      switch ifile
         case 1
            data.file1.w = w;
            data.file1.x = x;
         case 2
            data.file2.w = w;
            data.file2.x = x;
         otherwise
            disp('Error in read_data: cannot handle more than 2 files.')
            fprintf('nfile     = %d',nfile);
            fprintf('filenames = %s',filenames)
            return
      end
      fprintf('Write data to x%d and w%d\n',ifile,ifile)
      
      %read_transform_param
      
      %do_transform
      
      %if usereg
      fprintf('Finished reading %s\n',filelist(ifile).name);
   end
   fclose(fileID(ifile));
end

end