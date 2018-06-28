% Magnetic moments transformation from Cartesian to spherical Coord.
%
% Hongyang Zhou, hyzhou@umich.edu

%%
% Ask Xianzhe: dipole + induced magnetic field
% Here I try the way in Duling`s paper [2014]
% Table IV, Kivelson[2002]
Mz_dip = -727.3; My_dip = 52.8; Mx_dip = -18.4; %[nT]
% Table V, Kivelson[2002]
Mz = -716.8; My = 49.3; Mx = -22.2; %[nT]
phi = -7.1; %[degree]
A = 0.95; % amplitude factor

% G8
%Mx = -18.0; My = 51.8; Mz = -716.8;
%lambda = 285.9;

% G1
%Mx = -24.7; My = 82.5; Mz = -716.8;
lambda = 174.2;

% G2
%Mx = -29.3; My = 80.0; Mz = -716.8;
%lambda = 157.3;

% G7
%Mx = -20.9; My = 14.0; Mz = -716.8;
%lambda = 19.4;

% G28
%Mx = -19.3; My = 17.0; Mz = -716.8;
%lambda = 351.7;

% G29
%Mx= -18.4; My = 84.2; Mz = -716.8;

% G29
%B0x = -7; B0y = -78; B0z = -79; % [nT]
%lambda = 220.8; %[degree]


B0x = (-18)*sin((lambda-200)/180*pi);
B0y = (-86)*cos((lambda-200)/180*pi);

% There is a sign problem with the original eq. set!
Mx_ind = A*(-18)*sin((lambda+phi-200)/180*pi)*(0.5)
My_ind = A*(-86)*cos((lambda+phi-200)/180*pi)*(-0.5)
Mz_ind = 0;

Mx = Mx_dip + Mx_ind;
My = My_dip + My_ind;
Mz = Mz_dip + Mz_ind;

DipoleStrength = sqrt(Mx^2+My^2+Mz^2)

MagAxisThetaGeo = atan(sqrt(Mx^2+My^2)/abs(Mz))/pi*180

if Mz<0
   if My>0
      MagAxisPhiGeo = atan2(-My,-Mx)/pi*180+360
   else
      MagAxisPhiGeo = atan2(-My,-Mx)/pi*180
   end
else
   if My<0
      MagAxisPhiGeo = atan2(-My,-Mx)/pi*180+360
   else
      MagAxisPhiGeo = atan2(-My,-Mx)/pi*180
   end   
end
