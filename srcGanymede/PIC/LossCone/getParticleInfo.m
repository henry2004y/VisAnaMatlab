function [angle,Bx_P,By_P,Bz_P,B_P,particle,weight] = getParticleInfo
%GETPARTICLEINFO Calculate the derived particle information 
%
%
%OUTPUTS
% angle:    pitch angle for each particle, [nP,1]
% Bx_P,By_P,Bz_P: B field at particle locations, [nP,1]
% B_P:      B field strength at particle locations, [nP,1]
% particle: particle info inside the selected region [6,nP]
% weight:   particle weights,  [nP,1]

Species = Parameters.Species;

[xP,yP,zP,ux,uy,uz,weightP] = getParticle(Species);
[xF,yF,zF,Bx,By,Bz] = getField;

Region = Parameters.Region;
ncountmax = Parameters.ncountmax;
particle = Inf(6,ncountmax); weight = Inf(ncountmax,1);

nP = 0;
for iP=1:numel(xP)
   if xP(iP) >= Region(1) && xP(iP) <= Region(2) && ...
      yP(iP) >= Region(3) && yP(iP) <= Region(4) && ...
      zP(iP) >= Region(5) && zP(iP) <= Region(6)
      nP = nP + 1;
      particle(:,nP) = [xP(iP) yP(iP) zP(iP) ...
         ux(iP) uy(iP) uz(iP)]';
      weight(nP) = weightP(iP);
   end
end

% Cut-off the unused trails
particle = particle(:,1:nP); 
weight = weight(1:nP);

% Get pitch angle for each particle
Fx = griddedInterpolant(xF,yF,zF,Bx);
Fy = griddedInterpolant(xF,yF,zF,By);
Fz = griddedInterpolant(xF,yF,zF,Bz);

Bx_P = Fx(particle(1,:),particle(2,:),particle(3,:))';
By_P = Fy(particle(1,:),particle(2,:),particle(3,:))';
Bz_P = Fz(particle(1,:),particle(2,:),particle(3,:))';

B = [Bx_P By_P Bz_P];
U = particle(4:6,:)';
   
angle = atan2d(vecnorm(cross(B,U),2,2),dot(B,U,2));

% B Strength at particle positions
B_P = sqrt(Bx_P.^2 + By_P.^2 + Bz_P.^2);

% Plot pitch angles (weights included)
% Note: In order to see the "real" distribution, it is not correct to use
% angle directly as horizontal axis; instead, we should use cos(angle) as
% the true horizontal axis! Think of it about the circle.
[histw] = histwv(cosd(angle),abs(weight),-1,1,50);
figure
bar(linspace(0,180,50),histw)
% histogram(histw)

end
