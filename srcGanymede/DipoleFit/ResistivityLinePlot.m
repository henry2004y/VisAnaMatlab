% Resistivity line plot for Ganymede MHD-EPIC simulations
%
% PlotPub is used for better figure. Download it from Internet if you don't
% have it in your path.
%
% Hongyang Zhou, hyzhou@umich.edu 08/03/2017

clear; clc
%%
r   = [0.5 0.52 0.6 0.95 1.05 2];
%eta = [1 5e10 3e11 3e11 1 1]; % use 1 to replace 0 in log plot
eta = [1 6e9 6e11 6e11 1 1]; % use 1 to replace 0 in log plot
semilogy(r,eta,'*-','LineWidth',2); %axis tight
xlabel('radial distance [$R_G$]','Interpreter','LaTex');
ylabel('resistivity [$m^2/s$]','Interpreter','LaTex')
%set(gca,'FontSize',18)

%%
%opt.XLabel = 'radial distance r [R_G]'; % xlabel
%opt.YLabel = 'resistivity \eta [m^2/s]'; %ylabel

opt.Legend = {'Resistivity'};

setPlotProp(opt);