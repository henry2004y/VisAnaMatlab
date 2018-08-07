% Script for transforming SWMF PC outputs to vtk format
%
% Hongyang Zhou, hyzhou@umich.edu 03/24/2018

%% Pre-processing

flyby = 'G8';

switch flyby
   case 'G8'
      %filename = '3d_fluid_G8_350s.out';
      %filename = '~/Documents/research/Ganymede/data/3dPC_G8_t=1090.out';
      filename = '~/Documents/research/Ganymede/data/3d_fluid_G8_350s.out';
   case 'G28'
      %filename = '3d_fluid_G28_242s.out';
      filename = '~/Documents/research/Ganymede/data/3dPC_G28_t=580.out';
end

ipict = 1;

[filehead,data,list] = read_data(filename,'verbose',true,'npict',ipict);

x = data.file1.x(:,:,:,1);
y = data.file1.x(:,:,:,2);
z = data.file1.x(:,:,:,3);

% ndgrid --> meshgrid
x = permute(x,[2 1 3]);
y = permute(y,[2 1 3]);
z = permute(z,[2 1 3]);


func = 'Bx'; 
func_ = strcmpi(func,filehead.wnames);
Bx = data.file1.w(:,:,:,func_);
Bx = permute(Bx,[2 1 3]);

func = 'By'; 
func_ = strcmpi(func,filehead.wnames);
By = data.file1.w(:,:,:,func_);
By = permute(By,[2 1 3]);

func = 'Bz'; 
func_ = strcmpi(func,filehead.wnames);
Bz = data.file1.w(:,:,:,func_);
Bz = permute(Bz,[2 1 3]);

func = 'ps1'; 
func_ = strcmpi(func,filehead.wnames);
pi = data.file1.w(:,:,:,func_);
pi = permute(pi,[2 1 3]);


%% Save to vtk file
vtkwrite('G8_t350.vtk', 'structured_grid', x, y, z, ... 
   'vectors', 'magnetic_field', ...
   Bx, By, Bz, 'scalars', 'Pi', pi);




