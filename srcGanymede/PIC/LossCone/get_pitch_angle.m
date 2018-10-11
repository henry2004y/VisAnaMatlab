function [angle,B_P,particle,weight] = get_pitch_angle
%GET_PITCH_ANGLE
%
%INPUTS
%
%
%OUTPUTS
% angle: pitch angle for each particle
% B_P: surface B field strength for each particle
% particle: particle info inside the selected region
% weight: particle weights 

[xP,yP,zP,ux,uy,uz,weight] = get_particle('ion');
[xF,yF,zF,Bx,By,Bz] = get_field;

Region = Parameters.Region;
ncountmax = Parameters.ncountmax;
particle = Inf(6,ncountmax);

nP = 0;
for iP=1:numel(xP)
   if xP(iP) >= Region(1) && xP(iP) <= Region(2) && ...
      yP(iP) >= Region(3) && yP(iP) <= Region(4) && ...
      zP(iP) >= Region(5) && zP(iP) <= Region(6)
      %particle = [particle; ux(ipar) uy(ipar) uz(ipar)];
      nP = nP + 1;
      particle(:,nP) = [xP(iP) yP(iP) zP(iP) ...
         ux(iP) uy(iP) uz(iP)]';
   end
end

particle = particle(:,1:nP);
angle = Inf(nP,1);

% Get pitch angle for each particle 
Fx = griddedInterpolant(xF,yF,zF,Bx);
Fy = griddedInterpolant(xF,yF,zF,By);
Fz = griddedInterpolant(xF,yF,zF,Bz);

Bx_P = Inf(nP,1); By_P = Inf(nP,1); Bz_P = Inf(nP,1);
for iP=1:nP
   Bx_P(iP) = Fx(particle(1,iP),particle(2,iP),particle(3,iP));
   By_P(iP) = Fy(particle(1,iP),particle(2,iP),particle(3,iP));
   Bz_P(iP) = Fz(particle(1,iP),particle(2,iP),particle(3,iP));
   
   B = [Bx_P(iP) By_P(iP) Bz_P(iP)]; 
   U = [particle(4,iP) particle(5,iP) particle(6,iP)];
   
   angle(iP) = atan2d(norm(cross(B,U)),dot(B,U));   
end

% B Strength at particle positions
B_P = sqrt(Bx_P.^2 + By_P.^2 + Bz_P.^2);

% Plot
figure
histogram(angle)


end