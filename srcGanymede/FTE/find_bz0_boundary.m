function [ x3bc,y3bc,z3bc ] = find_bz0_boundary( filename,s,xThres )
%FIND_BZ0_BOUNDARY Returns the boundary points position for a volume.
%   
% Two conditions are used to pick the required boundary pts:
% 1. The distance from center should be larger than 1.6 Rg;
% 2. x must be negative (indicating upstream hemisphere).
%--------------------------------------------------------------------------

[head,data] = read_data(filename,'verbose',false);

data = data.file1;

x = data.x(:,:,:,1); 
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

VarIndex_ = strcmpi('bz',head.wnames);
bz = data.w(:,:,:,VarIndex_);

sig = sign(bz);
fx  = gradient(sig);
BCindex_ = find(fx~=0);

x3bc = x(BCindex_);
y3bc = y(BCindex_);
z3bc = z(BCindex_);

% Threshold for x
%xThres = -1.7;  % G8
%xThres = -1.75; % G28
% Find the outer boundary points
BCindex_ = x3bc < xThres;
x3bc = x3bc(BCindex_);
y3bc = y3bc(BCindex_);
z3bc = z3bc(BCindex_);



figure;scatter3(x3bc,y3bc,z3bc,'.'); axis equal

end