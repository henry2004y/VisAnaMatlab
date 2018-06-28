function [ centroid, boundaries, clusterNumber ] = find_FTE( var,varargin )
%FIND_FTE Find the number and boundaries of flux ropes
%% FTE counts
% # filter/smoothing
% # determine threshold for each variable
% # transform into binary(logical) matrix
% # Label blobs, measure blob properties
% # Merge close blobs, rule out small blobs
% # Identify moving blobs
%
% See references from image processing toolbox, and demos by experts.

%% Parameters and thresholds

p = inputParser;

defaultVarThreshold  = 7;
defaultAreaThreshold = 40;
defaultDistThreshold = 15;

addRequired(p,'var',@isnumeric)
addParameter(p,'VarThreshold',defaultVarThreshold,@isnumeric)
addParameter(p,'AreaThreshold',defaultAreaThreshold,@isnumeric)
addParameter(p,'DistThreshold',defaultDistThreshold,@isnumeric)

parse(p,var,varargin{:})

var = p.Results.var;
% Values larger than threshold is identified
VarThreshold = p.Results.VarThreshold;
% Blobs with area larger than threshold is identified
AreaThreshold = p.Results.AreaThreshold;
% Blobs with centroid distances less than threshold are merged into one
DistThreshold = p.Results.DistThreshold;


%% Build binary matrix
% I was thinking about using local maxima to identify peaks, but if the
% peak region does not have constant value, the result is not good.
% 
% Test for a new algorithm, Otsu`s method
% The input double type must be normalized to [0,1]. In this way I do not
% need to pick a threshold myself.
% VarThreshold = (VarThreshold-min(var(:))) ./ (max(var(:))-min(var(:)));
% var = (var-min(var(:))) ./ (max(var(:))-min(var(:)));
% VarThreshold = max(graythresh(var),VarThreshold);
% var_binary = imbinarize(var,VarThreshold);
%imshowpair(var,var_binary,'montage')

var_binary = logical(var);
var_binary(var <= VarThreshold) = false;
var_binary(var > VarThreshold)  = true;
var_binary(isnan(var)) = false;

% Using gradient info 
% var_grad = gradient(var);
% var_binary = logical(var);
% var_binary(var_grad <= 0.6) = false;
% var_binary(var_grad > 0.6) = true;

% Filtering out noise
% jblur = imgaussfilt(jv, 2);
% j_binary = jblur;
% j_binary(jblur<=threshold) = false;
% j_binary(jblur>threshold)  = true;
% j_binary(isnan(jblur)) = false;
% j_binary = logical(j_binary);

%%
% Label each blob so we can make measurements of it
% labeledImage is an integer-valued image where all pixels in the blobs 
% have values of 1, or 2, or 3, or ... etc.
blobMeasurements = regionprops(var_binary, 'centroid','PixelIdxList',...
   'PixelList','Area');

numberOfBlobs = size(blobMeasurements, 1);

% Rule out the small regions
for ib=1:numberOfBlobs
   if blobMeasurements(ib).Area < AreaThreshold
       var_binary(blobMeasurements(ib).PixelIdxList) = false;
   end
end

blobMeasurements = regionprops(var_binary, 'centroid','PixelIdxList',...
   'PixelList','Area');

boundaries = bwboundaries(var_binary,'noholes');
numberOfBoundaries = size(boundaries, 1);

% I need blob index and centroid info
% If the overlapped area is larger than a certain percent, treat them as 1
% blob.
% Maybe comparing the centroid distances is not a very good idea.

centroid = reshape([blobMeasurements.Centroid],[2,numberOfBoundaries])';
centriodDistance = pdist(centroid);

% Map that gives us which pixel belongs to which cluster
var_labeled = bwlabel(var_binary, 8);

% Create an array that tells us which object belongs to what cluster
membershipList = zeros(numberOfBoundaries,1);

% Initialize cluster count
clusterNumber = 1;

for ib = 1:numberOfBoundaries
   if membershipList(ib)==0
      membershipList(ib) = clusterNumber;  
      clusterNumber = clusterNumber + 1;
   end
   
   for jb = ib+1:numberOfBoundaries
      if centriodDistance((ib-1)*(numberOfBoundaries-ib/2)+jb-ib) ...
            <= DistThreshold        
         % Assign the same membership number
         membershipList(jb) = membershipList(ib);
      end
   end
   
end

% Give the actual number of blobs after merging
clusterNumber = clusterNumber - 1;

% Merge clusters by looping over the clusterNumber
for k = 1:clusterNumber
   blobIndex = find(membershipList==k);
   for ib = 1:length(blobIndex)
      placesToMark = blobMeasurements(blobIndex(ib)).PixelList;
      var_labeled(sub2ind(size(var),...
         placesToMark(:,2),placesToMark(:,1))) = k;
   end
end

end

