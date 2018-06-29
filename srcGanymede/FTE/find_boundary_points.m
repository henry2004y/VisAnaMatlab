function [ x3bc,y3bc,z3bc ] = find_boundary_points( filename,varargin )
%FIND_BOUNDARY_POINTS Returns the boundary points position for a volume.
%   
% Two conditions are used to pick the required boundary pts:
% # The distance from center should be larger than rThres;
% # x must be smaller than xThres (indicating upstream hemisphere).
%
% INPUT:
% filename: 3D output file name that contains status infomation
% s       : compact boundary factor (see built-in function *boundary*)
% DoPlot  : logical var, deciding doing scatter plot or not
% rThres  : radius threshold for picking boundary points
% xThres  : x coordinate threshold for picking boundary pts
%
% OUTPUT:
% x3bc,y3bc,z3bc: coordinates of boundary points
%
% Hongyang Zhou, hyzhou@umich.edu
% Modified 06/29/2018
%--------------------------------------------------------------------------

if nargin==0
   error('Not enough inputs.')
elseif nargin > 5
   error('Too Many Inputs.');
end

optargs = {0.5 true 1.5 -1.125}; % default parameters
optargs(1:nargin-1) = varargin;
% Place optional args in memorable variable names
[s, DoPlot, rThres, xThres] = optargs{:};


[head,data] = read_data(filename,'verbose',false);

data = data.file1;

x = data.x(:,:,:,1); 
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

VarIndex_ = strcmpi('status',head.wnames);
status = data.w(:,:,:,VarIndex_);

x3_1 = x(status==3 & x<0);
y3_1 = y(status==3 & x<0);
z3_1 = z(status==3 & x<0);

% Pick the compact boundary points coords.
k31  = boundary(x3_1,y3_1,z3_1,s);
bc31 = unique(k31);

x3bc = x3_1(bc31);
y3bc = y3_1(bc31);
z3bc = z3_1(bc31);

% Find the outer boundary points
mapindex_ = (x3bc.^2 + y3bc.^2 + z3bc.^2) > rThres^2 & x3bc < xThres;
x3bc = x3bc(mapindex_);
y3bc = y3bc(mapindex_);
z3bc = z3bc(mapindex_);

if DoPlot
   figure; scatter3(x3bc,y3bc,z3bc,'.'); axis equal
end

end

