function [ x3bc,y3bc,z3bc ] = find_bz0_boundary( filename,varargin )
%FIND_BZ0_BOUNDARY Returns the boundary points position for a volume.
%   
% Two conditions are used to pick the required boundary pts:
% 1. Bz changes sign;
% 2. x must be smaller than xThres (indicating upstream hemisphere).
%
% INPUT:
% filename: 3D output file name that contains status infomation
% DoPlot  : logical var, deciding doing scatter plot or not
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
elseif nargin > 3
   error('Too Many Inputs.');
end

optargs = {0.5 true -1.7}; % default parameters
% Threshold for x
%xThres = -1.7;  % G8
%xThres = -1.75; % G28
optargs(1:nargin-1) = varargin;
% Place optional args in memorable variable names
[DoPlot, xThres] = optargs{:};


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


% Find the outer boundary points
BCindex_ = x3bc < xThres;
x3bc = x3bc(BCindex_);
y3bc = y3bc(BCindex_);
z3bc = z3bc(BCindex_);


if DoPlot
   figure; scatter3(x3bc,y3bc,z3bc,'.'); axis equal
end
   
end