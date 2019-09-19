% Generalization of the Harris neutral current sheet configuration.
% The Coalescence of Magnetic Flux Ropes and Reconnecfion in the
% Magnetotail, Robert Richard, Raymond Walker, ...
% x-z plane
%
% Hongyang Zhou, hyzhou@umich.edu 09/02/2019

clear; close all; clc
%%
epsilon = 0.1; % controls the width of the islands and the sharpness of j
B0 = 1.0;    % magnetic field strength at large z

x = linspace(0,4*pi,100);
z = linspace(-1.2,1.2,100);

[X, Z] = meshgrid(x,z);

B0x = B0 .* sinh(Z)./(cosh(Z) + epsilon.*cos(X));
B0z = B0 .* epsilon.*sin(X)./(cosh(Z) + epsilon.*cos(X));

startx = linspace(0,4*pi,40);
startz = linspace(-1,1,40);
streamline(X, Z, B0x, B0z, startx, startz)
xlabel('x'); ylabel('z')
axis equal
set(gca, 'FontSize',14,'LineWidth',1.1)

