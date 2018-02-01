function [pictsize,filehead] = get_file_head( fileID,type )
%GET_FILE_HEAD Obtain the header information from BATSRUS output files.
% INPUT:
%   fileID: integer number allocated to file
%   type:   file type {'ascii', 'real4', 'binary', 'log'}
% OUTPUT:
%   pictsize: size of a single snapshot in bytes
%   filehead: file header info
%
%NOTE: An unformatted FORTRAN binary file contains record length 
% information along with each record. A record is a WRITE statement. That`s
% why I need to skip these records!
%--------------------------------------------------------------------------
persistent head

% Create a struct for substituting common block file_head
if isempty(head)
   % number of spatial dimensions
   field1 = 'ndim';       value1 = Inf;
   % 1st line often containing physcial unit names
   field2 = 'headline';   value2 = {''};
   % time step
   field3 = 'it';         value3 = Inf;
   % simulation time
   field4 = 'time';       value4 = Inf;
   % true for unstructured/non-Cartesian grids
   field5 = 'gencoord';   value5 = false;
   % number of scalar parameters
   field6 = 'neqpar';     value6 = int32(0);
   % number of variables
   field7 = 'nw';         value7 = Inf;
   % number of grid cells
   field8 = 'nx';         value8 = int32(0);
   % values of scalar parameters
   field9 = 'eqpar';      value9 = Inf;
   % array of coordinate/variable/param names
   field10 = 'variables'; value10 = {''};
   % array of variables names
   field11 = 'wnames';    value11 = {''};
   
   head = struct(field1,value1,field2,value2,field3,value3,...
      field4,value4,field5,value5,field6,value6,field7,value7,...
      field8,value8,field9,value9,field10,value10,field11, value11);
end

ftype = string(lower(type));

if strcmp(ftype,type); lenstr = 79; else lenstr=500; end

% Read header
pointer0 = ftell(fileID);

switch ftype
   case 'log'
      % This part is done in the main read_data function.
   case 'ascii'
      % This part needs to be tested! I think it is not usable now!
      head.headline = fgetl(fileID);
      line = fgetl(fileID);
      line = split(line);
      head.it = str2double(cell2mat(line(2)));
      head.time = str2double(cell2mat(line(3)));
      head.ndim = str2double(cell2mat(line(4)));
      head.neqpar = str2double(cell2mat(line(5)));
      head.nw = str2double(cell2mat(line(6)));
      head.gencoord = (head.ndim < 0);
      head.ndim = abs(head.ndim);
      line = str2double( split( fgetl(fileID) ) );      
      head.nx = line(2:end);
      if head.neqpar > 0
         line = str2double( split( fgetl(fileID) ) );
         head.eqpar = line(2:end);
      end
      varname = fgetl(fileID);
      
   case {'real4','binary'}
      fseek(fileID,4,'cof'); % skip record start tag.
      head.headline = fread(fileID,lenstr,'*char')';
      fseek(fileID,8,'cof'); % skip record end/start tags.
      head.it = fread(fileID,1,'*int');
      head.time = fread(fileID,1,'*float');
      head.ndim = fread(fileID,1,'*int');
      head.gencoord = (head.ndim < 0);
      head.ndim = abs(head.ndim);
      head.neqpar = fread(fileID,1,'*int');
      head.nw = fread(fileID,1,'*int');
      fseek(fileID,8,'cof'); % skip record end/start tags.
      head.nx = fread(fileID,head.ndim,'*int');
      fseek(fileID,8,'cof'); % skip record end/start tags.
      if head.neqpar > 0
         head.eqpar = fread(fileID,head.neqpar,'*float');
      end
      fseek(fileID,8,'cof'); % skip record end/start tags.
      varname = fread(fileID,lenstr,'*char')';
      fseek(fileID,4,'cof'); % skip record end tag.
end

%Header length
pointer1 = ftell(fileID);
headlen = pointer1-pointer0;

% Calculate the snapshot size = header + data + recordmarks
nxs = prod(head.nx);

switch ftype
   case 'log'
      pictsize = 1;
   case 'ascii'
      pictsize = headlen + (18*(head.ndim+head.nw)+1)*nxs;
   case 'binary'
      pictsize = headlen + 8*(1+head.nw)+8*(head.ndim+head.nw)*nxs;
   case 'real4'
      pictsize = headlen + 8*(1+head.nw) + 4*(head.ndim+head.nw)*nxs;
end

if nargout>1
   % Set variables array
   %head.variables = strsplit(varname); % returns a cell array
   head.variables = split(varname);     % returns a string array
      
   if ftype == 'real4' | ftype == 'binary'
      % remove the last null string
      head.variables = head.variables(1:end-1);
   end
   % Return head as filehead
   filehead = head;
end

end