function [dBx,dBy,dBz,limits] = GetMeanField(Dir,fnameParticle,...
   fnameField,limits)
%GetMeanField Get the average field direction in limited region
%   * Obtain the limits from particle data
%   * Extract the average field from field data

%% Obtain the limits
if nargin < 4
   filename = fullfile(Dir,fnameParticle);
   [filehead,data] = read_data(filename,'verbose',false);
   
   data = data.file1;
   
   xmin = min(data.x(:,:,:,1));
   xmax = max(data.x(:,:,:,1));
   ymin = min(data.x(:,:,:,2));
   ymax = max(data.x(:,:,:,2));
   zmin = min(data.x(:,:,:,3));
   zmax = max(data.x(:,:,:,3));
   
   % xmin = -1.9;
   % xmax = -1.85;
   % ymin = -0.05;
   % ymax = 0.05;
   % zmin = -0.60;
   % zmax = -0.55;
   
   if strcmpi(filehead.headline(1:2),'SI')
      limits = [xmin,xmax,ymin,ymax,zmin,zmax]/2634000;
   else % planetary units
      limits = [xmin,xmax,ymin,ymax,zmin,zmax];
   end
end   
   
%% Get the average field direction in limited region
filename= fullfile(Dir,fnameField);
[filehead,data] = read_data(filename,'verbose',false);

data = data.file1;

x = data.x(:,:,:,1);
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

x = permute(x,[2 1 3]);
y = permute(y,[2 1 3]);
z = permute(z,[2 1 3]);

Bx_ = strcmpi('Bx',filehead.wnames);
Bx = data.w(:,:,:,Bx_);
Bx = permute(Bx,[2 1 3]);

By_ = strcmpi('By',filehead.wnames);
By = data.w(:,:,:,By_);
By = permute(By,[2 1 3]);

Bz_ = strcmpi('Bz',filehead.wnames);
Bz = data.w(:,:,:,Bz_);
Bz = permute(Bz,[2 1 3]);

[~,~,~,Bx,By,Bz] = subvolume(x,y,z,Bx,By,Bz,limits);

% Average over the selected volume
Bx = mean(Bx(:)); By = mean(By(:)); Bz = mean(Bz(:));

% Unify vector
Length = norm([Bx By Bz]);
dBx = Bx/Length; dBy = By/Length; dBz = Bz/Length;

end

