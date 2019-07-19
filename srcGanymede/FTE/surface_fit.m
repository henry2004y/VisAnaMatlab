function [fitresult,gof] = surface_fit(x3bc,y3bc,z3bc,varargin)
%SURFACE_FIT Fit the closed field line boundary with hypersurface fit
%
% INPUTS:
% x3bc,y3bc,z3bc: 1D array of coordinates
% TypeFit: fit model type (see MATLAB document)
% DoPlot : logical, display figure or not
%
% OUTPUTS:
% fitresult: fit model
% gof: goodess of fit
%
% Hongyang Zhou, hyzhou@umich.edu 06/29/2018

if nargin<3
   error('Not enough input arguments.')
end

optargs = {'poly55' false}; % default parameters
optargs(1:nargin-3) = varargin;
[TypeFit,DoPlot] = optargs{:};

% Set up fittype and options.
ft = fittype( TypeFit );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

if DoPlot
   % Plot fit with data.
   figure(1); %hold on
   h = plot( fitresult );
   legend(h, 'poly5, x=x(y,z)', 'Location', 'NorthEast');
   % Label axes
   xlabel('x [R_G]')
   ylabel('y [R_G]')
   zlabel('z [R_G]')
   grid on; axis equal
   
   xx = get(h, 'XData');
   yy = get(h, 'YData');
   zz = get(h, 'Zdata');
   set(h, 'XData', zz, 'YData', xx, 'ZData', yy);
   
   hold on;
   scatter3(x3bc,y3bc,z3bc,20,'r','filled'); hold off
   %axis tight
   xlim([-2 0]); ylim([-1.5 1.5]); zlim([-0.6 0.8]);
   set(gca,'FontSize',14,'LineWidth',1.2)
end

end
