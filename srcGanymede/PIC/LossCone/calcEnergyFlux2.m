function CalcEnergyFlux2(particle,angle,phi1,theta1,Bsurf,B_P)
%CALCENERGYFLUX2 Calc energetic fluxes at the surface
% Calculate particle velocities on the surface, and then compute the
% moments and fluxes.
% The problem is: how can I get the density at the surface?
%
%INPUTS:
% particle: positions, velocities and weights, [7,nP]
% angle: pitch angle
% phi1:
% theta1:
% Bsurf:
% B_P
%

if strcmp(Parameters.Species,'ion')
   m = Parameters.mi * Parameters.mp;
else
   m = Parameters.me;
end

% Perpendicular velocities at the current location
vPerp = sqrt(sum(particle(4:6,:).^2,1)') .* sind(angle);
% Mapping particles to the surface
vPerp = sqrt(Bsurf./B_P) .* vPerp;
vPar  = sqrt(sum(particle(4:6,:).^2,1)' - vPerp.^2);

% Calculate average velocities by bins
bins = 30;
phiMin = min(phi1); phiMax = max(phi1);
thetaMin = min(theta1); thetaMax = max(theta1);
dPhi = (phiMax - phiMin)/(bins - 1); 
dTheta = (thetaMax - thetaMin)/(bins - 1);
subs = [round((phi1-phiMin)/dPhi)+1 round((theta1-thetaMin)/dTheta)+1]; 

% uPar = accumarray(subs,vPar.*particle(7,:)') ./ ...
%    accumarray(subs,particle(7,:)');

% Generate mesh
theta = thetaMin:dTheta:thetaMax;
phi = phiMin:dPhi:phiMax;
[Phi,Theta] = meshgrid(phi,theta);

flux = accumarray(subs,(0.5*m*(vPar.^2 + vPerp.^2).*vPar.*particle(7,:)'));

% Get the flux normal to the surface
[FBxSurf,FBySurf,FBzSurf] = getBsurface(false);
BxSurf = FBxSurf(Phi',Theta');
BySurf = FBySurf(Phi',Theta');
BzSurf = FBzSurf(Phi',Theta');

Br = BxSurf' .* cosd(Phi) .* cosd(Theta) + ...
     BySurf' .* sind(Phi) .* cosd(Theta) + ...
     BzSurf' .* sind(Theta);
  
flux = flux .* abs(Br) ./ sqrt(BxSurf.^2 + BySurf.^2 + BzSurf.^2);

%
figure
axesm('ortho','origin',[45 45]); 
axis off;
gridm on; 
framem on;
mlabel('equator')
setm(gca,'Origin',[0 180 0])
plabel(120); 
plabel('fontweight','bold')

h1 = surfm(Theta,Phi,flux'); %colorbar

%
% figure
% axesm('ortho','origin',[45 45]); 
% axis off;
% gridm on; 
% framem on;
% mlabel('equator')
% plabel(0); 
% plabel('fontweight','bold')
% 
% [n,c] = hist3([phi1,theta1],'CdataMode','auto','Nbins',[30 30]);
% h1 = surfm(Theta,Phi,n');
% setm(gca,'Origin',[0 180 0])

end

