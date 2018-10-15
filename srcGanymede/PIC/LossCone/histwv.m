function [histw, histv] = histwv(v, w, Vmin, Vmax, bins) 
%HISTWV Generate data for weighted histogram
%
%Inputs: 
% v     - values 
% w     - weights 
% Vmin  - minimum value 
% Vmax  - max value 
% bins  - number of bins (inclusive) 
%
%Outputs: 
% histw - weighted histogram 
% histv (optional) - histogram of values 

delta = (Vmax - Vmin)/(bins - 1); 
subs = round((v - Vmin)/delta)+1; 

histv = accumarray(subs,1,[bins,1]); 
histw = accumarray(subs,w,[bins,1]); 
end