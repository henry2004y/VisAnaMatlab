function [ x,w ] = get_pict_bin( fileID,filehead )
%get_pict_bin Read binary format data
%   
%--------------------------------------------------------------------------

% Read coordinates and values
switch filehead.ndim
   case 1 % 1D
      n1 = filehead.nx(1);
      fseek(fileID,4,'cof'); % skip record start tag.
      x  = fread(fileID,[n1,filehead.ndim],'*double');
      w  = Inf(n1,filehead.nw);
      fseek(fileID,8,'cof'); % skip record end/start tags.
      for iw=1:filehead.nw
         w(:,iw) =fread(fileID,n1,'*double');
         fseek(fileID,8,'cof'); % skip record end/start tags.
      end
   case 2 % 2D
      n1 = filehead.nx(1);
      n2 = filehead.nx(2);
      fseek(fileID,4,'cof'); % skip record start tag.
      x  = fread(fileID,n1*n2*filehead.ndim,'*double');
      x  = reshape(x,[n1,n2,filehead.ndim]);
      w  = Inf(n1,n2,filehead.nw);
      fseek(fileID,8,'cof'); % skip record end/start tags.
      for iw=1:filehead.nw
         w(:,:,iw) =fread(fileID,[n1,n2],'*double');
         fseek(fileID,8,'cof'); % skip record end/start tags.
      end      
   case 3 % 3D
      n1 = filehead.nx(1);
      n2 = filehead.nx(2);
      n3 = filehead.nx(3);
      fseek(fileID,4,'cof'); % skip record start tag.
      x  = fread(fileID,n1*n2*n3*filehead.ndim,'*double');
      w  = Inf(n1,n2,n3,filehead.nw);
      fseek(fileID,8,'cof'); % skip record end/start tags.
      for iw=1:filehead.nw
         w(:,:,:,iw) =fread(fileID,[n1,n2,n3],'*double');
         fseek(fileID,8,'cof'); % skip record end/start tags.
      end
end

end

