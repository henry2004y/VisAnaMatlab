function [Fx,Fy,Fz] = get_Bsurface(DoPlot)
%GET_BSURFACE
%
%INPUTS
% DoPlot: Plot the surface B field
%
%OUTPUTS
% Fx,Fy,Fz: griddedinterpolant of B field on the 2D surface


Dir = '~/Ganymede/MOP2018/runG8_PIC_1200s/EnergeticFlux';
fname = 'shl_var_1_t00000557_n00250489.out';

[filehead,data] = read_data(fullfile(Dir,fname));
data = data.file1;

r = data.x(2,1,1,1); % cut radius [Rg]
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

Bx = squeeze(Bx(2,:,:));
By = squeeze(By(2,:,:));
Bz = squeeze(Bz(2,:,:));

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