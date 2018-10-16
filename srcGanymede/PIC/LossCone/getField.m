function [xF,yF,zF,Bx,By,Bz] = getField
%GETFIELD Read B field info from PIC output
%   
%INPUT
%
%OUTPUT:
% xF,yF,zF: grid locations
% Bx,By,Bz: B fields

% Set parameters
Dir = Parameters.Dir;
fnameField = Parameters.fnameField;

% Field data
[filehead,data] = read_data(fullfile(Dir,fnameField),'verbose',false);

data = data.file1;

bx_ = strcmpi('bx',filehead.wnames);
by_ = strcmpi('by',filehead.wnames);
bz_ = strcmpi('bz',filehead.wnames);

xF = data.x(:,:,:,1);       % [Rg]
yF = data.x(:,:,:,2);       % [Rg]
zF = data.x(:,:,:,3);       % [Rg]
Bx = data.w(:,:,:,bx_);      % [nT]
By = data.w(:,:,:,by_);
Bz = data.w(:,:,:,bz_);

end

