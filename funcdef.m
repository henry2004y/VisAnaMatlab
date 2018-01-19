function [ output_args ] = funcdef( input_args )
%FUNCDEF Generate derived variables from original variables for plotting.
%   Ideas borrowed from IDL version
%   This is a function called by plot_data and animate_data.
%   xx:    array contains the 'ndim' components of the coordinates.
%   w:     array contains the 'nw' variables.
%   func:  string describes the function to be returned.
%
%   Expressions formed from
%   1. standard variables:   rho,ux,uy,uz,p,bx,by,bz
%   2. standard coordinates: x,y,z,r
%   3. scalar parameters:    gamma, gamme,rbody,c0,mi,me
%
%--------------------------------------------------------------------------

% Define various functions of the basic MHD variables
% The function names are evaluated in lower case

% Variable names
wnames = strlowecse(variables());

% Set the coordinate arrays x,y,z if they occur in the variable names


% Set the variables from w
for iw=0:nw-1
   switch ndim
      case 1 
         switch wnames(iw) 
            case 'rho' 
               rho = w(:,iw);
               
         end
      case 2
         switch wnames(is)
            
         end
      case 3
         
   end
end

% Extra variables
uu = ux^2 + uy^2 + uz^2;
bb = bx^2 + by^2 + bz^2;
u  = sqrt(uu);           % speed
b  = sqrt(bb);           % magnetic field strength

% Calculate x1,x2,...,z3 field aligned basis vectors if needed


% Add functions to the basic variable list


end

