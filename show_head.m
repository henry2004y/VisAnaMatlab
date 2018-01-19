function show_head(file,ifile,filehead)
%show_head Displaying header information of ifile
%   Detailed explanation goes here
%--------------------------------------------------------------------------   
   
   disp('----------------------')
   fprintf('ifile     = %d\n',ifile)
   fprintf('filename  = %s\n',file.name)
   fprintf('filetype  = %s\n',file.type)
   fprintf('headline  = %s\n',filehead.headline)
   fprintf('it        = %d\n',filehead.it)
   fprintf('time      = %f\n',filehead.time)
   fprintf('gencoord  = %d\n',filehead.gencoord)
   fprintf('ndim      = %d\n',filehead.ndim)
   fprintf('neqpar    = %d\n',filehead.neqpar)
   fprintf('nw        = %d\n',filehead.nw)
   fprintf('nx        = %d\n',filehead.nx)
   disp('----------------------')
   
   if filehead.neqpar > 0
      fprintf('parameters = %15.8f\n',filehead.eqpar)
      fprintf('coord names= %s\n',filehead.variables{1:filehead.ndim})
      fprintf('var   names= %s\n',...
      filehead.variables{filehead.ndim+1:filehead.ndim+filehead.nw})
      fprintf('param names= %s\n',...
         filehead.variables{filehead.ndim+filehead.nw+1:end})
      disp('=======================')
   end
   
end