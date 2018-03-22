function [ filehead,data ] = read_log_data( filename )
%read_log_data Read information from log file.
%
%--------------------------------------------------------------------------

data = importdata(filename);

if isstruct(data)
   if isfield(data,'colheaders')
      filehead.headline  = data.textdata{1};
      filehead.variables = data.colheaders;
      filehead.ndim      = 1;
      filehead.it        = 0;
      filehead.time      = 0.0;
      filehead.gencoord  = false;
      filehead.nx        = 1;
      filehead.nw        = numel(data.colheaders);
   else
      filehead.variables = split(data.textdata);
   end
   
   data = data.data;
else
   filehead = [];
end


end