function p = streamlinec(s, X, Y, Z, C)
%STREAMLINEC Plot the streamlines with colors specified along the lines.
%INPUT:
% s: original streamline function
% X: grid points X
% Y: grid points Y
% Z: grid points Z
%OUTPUT:
% p: line handles, which will be a numeric value. You can control the
% properties of patches using p = findobj(p).
%
% Author: Matthew Crema, modified by Hongyang Zhou, 08/20/2017
%--------------------------------------------------------------------------

% Allocate space for output array of handles
p = Inf(length(s), 1);

% Iterate over streamlines
for stline = 1:length(s)
   % Get all vertices and faces for the current streamline
   v = s{stline};
   % This is amazing: faces can be lines!
   f = [(1:length(v)-1)' , (2:length(v))'];
   
   % Interpolate the color field at each vertex of the current streamline
   fvcdata = real(interp3(X, Y, Z, C, v(:,1), v(:,2), v(:,3)));
   
   p(stline) = patch('vertices', v, 'faces', f, ...
      'facevertexcdata', fvcdata, 'edgecolor', 'interp', ...
      'facecolor', 'none', 'facelighting', 'phong', ...
      'visible','on');
end

%set(p,'visible','on')
% Return handles in a cell array (hyzhou: which is not needed)
%p = {p};

end

