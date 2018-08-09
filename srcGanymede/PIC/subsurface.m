function [newx, newy, newdata] = subsurface(varargin)
%SUBSURFACE Extract subset of surface dataset.
%  This is a simplified version of subvolume.

x      = varargin{1};
y      = varargin{2};
data   = varargin{3};
limits = varargin{4};

if numel(limits)~=4
  error('Reduction must be [xmin xmax ymin ymax]');
end

if limits(1) > limits(2)
  error(message('MATLAB:subvolume:InvalidReductionXRange'));
end
if limits(3) > limits(4)
  error(message('MATLAB:subvolume:InvalidReductionYRange'));
end

sz = size(data);

hx = x(:,1);
hy = y(1,:);

if isnan(limits(1)),  limits(1) = min(hx); end
if isnan(limits(3)),  limits(3) = min(hy); end
if isnan(limits(2)),  limits(2) = max(hx); end
if isnan(limits(4)),  limits(4) = max(hy); end

xind = find(limits(1)<=hx & hx<=limits(2));
yind = find(limits(3)<=hy & hy<=limits(4));

newdata = subdata(data, xind, yind, sz);

newx = x(xind, yind);
newy = y(xind, yind);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newdata = subdata(data, xind, yind, sz)
newdata = data(xind, yind);
newsz = size(newdata);
if length(sz)>2
  newdata = reshape(newdata, [newsz(1:3) sz(4:end)]);
end

end

end

