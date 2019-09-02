function [ x,w ] = get_pict_asc( fileID,filehead )
%get_pict_asc Read ascii format data
%   
%--------------------------------------------------------------------------

ndim = filehead.ndim;
nw   = filehead.nw;

% Read coordinates and values row by row
switch ndim
   case 1 % 1D
      n1 = filehead.nx(1); 
      x  = Inf(n1,ndim);
      w  = Inf(n1,nw);
      for ix=1:n1
         temp = str2double(split(fgetl(fileID)));
         x(ix,:) = temp(2);
         w(ix,:) = temp(3:end);
      end
   case 2 % 2D
      n1 = filehead.nx(1);
      n2 = filehead.nx(2);
      x  = Inf(n1,n2,ndim);
      w  = Inf(n1,n2,nw);
      for ix1=1:n1
         for ix2=1:n2
            temp = str2double(split(fgetl(fileID)));
            x(ix1,ix2,:) = temp(2);
            w(ix1,ix2,:) = temp(3:end);
         end
      end      
   case 3 % 3D
      % This triple for loop is very slow.
      n1 = filehead.nx(1);
      n2 = filehead.nx(2);
      n3 = filehead.nx(3);
      x  = Inf(n1,n2,n3,ndim);
      w  = Inf(n1,n2,n3,nw);
      for ix1=1:n1
         for ix2=1:n2
            for ix3=1:n3
               temp = str2double(split(fgetl(fileID)));
               x(ix1,ix2,ix3,:) = temp(2);
               w(ix1,ix2,ix3,:) = temp(3:end);
            end
         end
      end
      
end

end
