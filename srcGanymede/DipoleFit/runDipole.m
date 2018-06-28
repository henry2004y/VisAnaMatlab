% Calculate dipole moments for Ganymede`s simulation input
% Hongyang Zhou, hyzhou@umich.edu 05/03/2017
%---------------------------------------------------------

% Note: I realized one problem: for Ganymede, we use GPhiO coordinates,
% while in BATSRUS we use GSE coordinates. The x and y axis are different!
% Maybe this is the problem for our simulation!
% 1. x,y --> -x,-y
% 2. x,y --> y,-x

%% Dipole fit from Kivelson [2002], Table V
Mx0 = -22.2; % [nT]
My0 = 49.3;  % [nT]
%Mx0 = 49.3; % [nT]
%My0 = 22.2;  % [nT]
Mz0 = -716.8;% [nT]
alpha = 0.84; % ratio to a perfect conducting sphere

% Background magnetic field for each flyby
%G1
Bbk1 = [6,-79,-79];  % [nT]
%G2
Bbk2 = [17,-73,-85]; % [nT]
%G7
Bbk7 = [-3,84,-76];  % [nT]
%G8
%Bbk8 = [-11,11,-77]; % [nT]
%Bbk8 = [-11,-6,-77]; % [nT]
Bbk8 = [-10,-6,-86]; % [nT]
%G28
Bbk28 = [-7,78,-76]; % [nT]
%G29
Bbk29 = [-9,-83,-79]; % [nT]

% Magnetic moment for each flyby as a sum of permanent dipole and 
% induced dipole from background magnetic field
%G1
Mx1 = Mx0-alpha*0.5*Bbk1(1);
My1 = My0-alpha*0.5*Bbk1(2);
Mz1 = Mz0;

%G2
Mx2 = Mx0-alpha*0.5*Bbk2(1);
My2 = My0-alpha*0.5*Bbk2(2);
Mz2 = Mz0;

%G7
Mx7 = Mx0-alpha*0.5*Bbk7(1);
My7 = My0-alpha*0.5*Bbk7(2);
Mz7 = Mz0;

%G8
Mx8 = Mx0-alpha*0.5*Bbk8(1);
My8 = My0-alpha*0.5*Bbk8(2);
Mz8 = Mz0;

%G28
Mx28 = Mx0-alpha*0.5*Bbk28(1);
My28 = My0-alpha*0.5*Bbk28(2);
Mz28 = Mz0;

%G29
Mx29 = Mx0-alpha*0.5*Bbk29(1);
My29 = My0-alpha*0.5*Bbk29(2);
Mz29 = Mz0;


% Transform from Cartesian to Spherical Coord.
%G1
[Theta1,Phi1,S1] = dipole_cart2sph(Mx1,My1,Mz1);
%G2
[Theta2,Phi2,S2] = dipole_cart2sph(Mx2,My2,Mz2);
%G7
[Theta7,Phi7,S7] = dipole_cart2sph(Mx7,My7,Mz7);
%G8
[Theta8,Phi8,S8] = dipole_cart2sph(Mx8,My8,Mz8);
%G28
[Theta28,Phi28,S28] = dipole_cart2sph(Mx28,My28,Mz29);
%G29
[Theta29,Phi29,S29] = dipole_cart2sph(Mx29,My29,Mz29);

%% Create a table for viewing
Phi = [Phi1,Phi2,Phi7,Phi8,Phi28,Phi29];
Theta = [Theta1,Theta2,Theta7,Theta8,Theta28,Theta29];
S = [S1,S2,S7,S8,S28,S29];
flyby = {'G1','G2','G7','G8','G28','G29'};
title = {'Theta','Phi','Strength'};

T = table(Theta',Phi',S','VariableNames',title,'RowNames',flyby)

% Save table into Excel format
%writetable(T,'DipoleAngle','FileType','spreadsheet')