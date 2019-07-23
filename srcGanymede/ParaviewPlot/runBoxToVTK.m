% Script for transforming SWMF GM box outputs to vtk format
%
% Hongyang Zhou, hyzhou@umich.edu 07/20/2019

clear; clc
%% Pre-processing

filenames = dir('box*');

for ifile = 1:length(filenames)
   
   filename = filenames(ifile).name;
   
   ipict = 1;
   
   [filehead,data,list] = read_data(filename,'verbose',true,'npict',ipict);
   
   data = data.file1;
   
   x = data.x(:,:,:,1);
   y = data.x(:,:,:,2);
   z = data.x(:,:,:,3);
   
   % ndgrid --> meshgrid
   x = permute(x,[2 1 3]);
   y = permute(y,[2 1 3]);
   z = permute(z,[2 1 3]);
   
   
   func = 'Bx';
   func_ = strcmpi(func,filehead.wnames);
   Bx = data.w(:,:,:,func_);
   Bx = permute(Bx,[2 1 3]);
   
   func = 'By';
   func_ = strcmpi(func,filehead.wnames);
   By = data.w(:,:,:,func_);
   By = permute(By,[2 1 3]);
   
   func = 'Bz';
   func_ = strcmpi(func,filehead.wnames);
   Bz = data.w(:,:,:,func_);
   Bz = permute(Bz,[2 1 3]);
   
   func = 'p';
   func_ = strcmpi(func,filehead.wnames);
   p = data.w(:,:,:,func_);
   p = permute(p,[2 1 3]);
   
   % Save to vtk file
   outname = strcat(filename(1:19),'.vtk');
   
   vtkwrite(outname, 'structured_grid', x, y, z, ...
      'vectors', 'B', Bx, By, Bz, 'scalars', 'P', p, 'binary');
end
