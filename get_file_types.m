function [filelist,fileID,pictsize] = get_file_types(nfile,filelist)
%get_file_types Get the type of output files
%  
%INPUTS:
% nfile: number of files
% filelist: struct of files that need to be fulfilled.
%OUTPUTS:
% filelist: completed file structs
% fileID: file # for accessing data
% pictsize: size (in bytes) of one snapshot
%--------------------------------------------------------------------------

fileID   = Inf(nfile,1);
pictsize = Inf(nfile,1);

for ifile=1:nfile
   f = fullfile(filelist(ifile).folder,filelist(ifile).name);
   fileID(ifile) = fopen(f,'r');
   
   % Check the appendix of file names
   % I realized that this may not be a robust way, especially in linux,
   % because you can have any appendix you want!!! Gabor uses a trick: the 
   % first 4 bytes decides the file type
   switch string(filelist(ifile).name(...
         strfind(filelist(ifile).name,'.')+1 : end))
      case {'log','sat'}
         filelist(ifile).type   = 'log';
         filelist(ifile).npictinfiles = 1;
                 
      otherwise
         % Obtain filetype based on the length info in the first 4 bytes
         lenhead = fread(fileID(ifile),1,'long');
         
         if lenhead~=79 && lenhead~=500
            filelist(ifile).type = 'ascii';
         else
            % The length of the 2nd line decides between real4 and binary
            % since it contains the time, which is real*8 or real*4
            head = fread(fileID(ifile),lenhead+4,'*char')';
            len = fread(fileID(ifile),1,'long');
            switch len
               case 20
                  filelist(ifile).type = 'real4';
               case 24
                  filelist(ifile).type = 'binary';
               otherwise
                  error(['Error in get_file_types: '...
                     'strange unformatted file:%s'],filelist(ifile).name)
            end
            if lenhead==500
               filelist(ifile).type = upper(filelist(ifile).type); 
            end
         end

      % Obtain file size and number of snapshots
      frewind(fileID(ifile));
      [pictsize(ifile)] = ...
         get_file_head(fileID(ifile),filelist(ifile).type);
      filelist(ifile).npictinfiles = ...
         floor(filelist(ifile).bytes / pictsize(ifile));
   end 
end

end