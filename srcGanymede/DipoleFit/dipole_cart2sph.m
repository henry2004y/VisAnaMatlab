function [ MagAxisThetaGeo, MagAxisPhiGeo, DipoleStrength ] = ...
   dipole_cart2sph( Mx,My,Mz )
% Magnetic moments transformation from Cartesian to spherical Coord.
%
%INPUTS:
% Mx: -h11 [nT]
% My:  g11 [nT]
% Mz:  g01 [nT]
%
%OUTPUTS:
% MaxAxisThetaGeo: angle from z axis in GSM [degree]
% MaxAxisPhiGeo  : angle from x axis in GSM [degree]
% DipoleStrength : dipole strength [nT]
%
% Hongyang Zhou, hyzhou@umich.edu
%-------------------------------------------------------------------

DipoleStrength = sqrt(Mx^2+My^2+Mz^2);

MagAxisThetaGeo = atan(sqrt(Mx^2+My^2)/abs(Mz))/pi*180;

if Mz<0
   if My>0
      MagAxisPhiGeo = atan2(-My,-Mx)/pi*180+360;
   else
      MagAxisPhiGeo = atan2(-My,-Mx)/pi*180;
   end
else
   if My<0
      MagAxisPhiGeo = atan2(-My,-Mx)/pi*180+360;
   else
      MagAxisPhiGeo = atan2(-My,-Mx)/pi*180;
   end   
end

end

