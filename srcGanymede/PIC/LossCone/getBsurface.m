function [Fx,Fy,Fz] = getBsurface(DoPlot)
%GETBSURFACE
%
%INPUTS
% DoPlot: Plot the surface B field
%
%OUTPUTS
% Fx,Fy,Fz: griddedinterpolant of B field on the 2D surface

% Set parameters
Dir = Parameters.Dir;
fname = Parameters.fnameSurf;


[filehead,data] = read_data(fullfile(Dir,fname),'verbose',false);
data = data.file1;

rIndex_ = 1;

r = data.x(rIndex_,1,1,1); % cut radius [Rg]
lon = data.x(:,:,:,2);
lat = data.x(:,:,:,3);

lat = squeeze(lat(2,:,:)); 
lon = squeeze(lon(2,:,:)); 

bx_ = strcmpi('bx',filehead.wnames);
by_ = strcmpi('by',filehead.wnames);
bz_ = strcmpi('bz',filehead.wnames);

Bx = data.w(:,:,:,bx_);      % [nT]
By = data.w(:,:,:,by_);
Bz = data.w(:,:,:,bz_);

Bx = squeeze(Bx(rIndex_,:,:));
By = squeeze(By(rIndex_,:,:));
Bz = squeeze(Bz(rIndex_,:,:));

Fx = griddedInterpolant(lon,lat,Bx);
Fy = griddedInterpolant(lon,lat,By);
Fz = griddedInterpolant(lon,lat,Bz);

if nargin > 0 && DoPlot
   % ndgrid --> meshgrid for plotting
   lat = permute(lat,[2 1]);
   lon = permute(lon,[2 1]);
   Bx = permute(Bx,[2 1]);
   By = permute(By,[2 1]);
   Bz = permute(Bz,[2 1]);
   
   figure
   subplot(221)
   contourf(lon,lat,Bx)
   title('Bx [nT]'); colorbar
   ylabel('Latitude')
   subplot(222)
   contourf(lon,lat,By)
   title('By [nT]'); colorbar
   subplot(223)
   contourf(lon,lat,Bz)
   title('Bz [nT]'); colorbar
   xlabel('Longitude'); ylabel('Latitude')
   subplot(224)
   contourf(lon,lat,sqrt(Bx.^2+By.^2+Bz.^2))
   title('B [nT]'); colorbar
   xlabel('Longitude')
end

end
