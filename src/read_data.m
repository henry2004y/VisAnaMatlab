function [ filehead,data,filelist ] = read_data( filename, varargin )
%read_data Read data from BATSRUS output files
%   Read the npict-th snapshot from an ascii or binary data file into          
%   the x (coordinates) and w (data) arrays.
%
%   filehead = read_data(filename)
%
%   [filehead, data] = read_data(filename,'npict',10)
%
%   [filehead, data, filelist] = read_data(filename, 'npict', 2, 
%      'verbose', false)
%                                                                   
%
%   The "x" and "w" arrays and the header info will be read from the file.
%
%   The same npict-th snapshot can be read from up to 10 files by e.g.
%   setting                                                                       
% filename='data/file1.ini data/file2.out'
%
%   In this case the data is read into x0,w0 and x1,w1 for the two files
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
