function getTemperature(particle,uBulk)
%GETTEMPERATURE Calculate the energetic fluxes at the surface
% First calculate the energetic fluxes along the field line at the original
% locations, and then mapping to the surface.
%
%INPUTS:
% particle: particle positions, velocities and weights, [7,nP]
%
%OUTPUT:
% Te: electron temperature in [eV]



%% Calculate at the original locations

q         = Parameters.q;
kB        = Parameters.kB;
me        = Parameters.me;
mp        = Parameters.mp;
mi        = Parameters.mi;
No2SiMass = Parameters.No2SiMass;

Te = 0.5*sum(sum(particle(4:6,:).^2,1).*particle(7,:))*No2SiMass / ...
   (size(particle,2)*0.03*No2SiMass/me) / q

Te = 0.5*sum(sum(particle(4:6,:).^2,1).*particle(7,:))*No2SiMass / ...
   (size(particle,2)*No2SiMass/me) / kB

Te1 = 0.5*sum(sum((particle(4:6,:) - uBulk).^2,1)...
   .*particle(7,:))*No2SiMass / ...
   (size(particle,2)*0.5*No2SiMass/me) / q

Ethres = 14.3; %[eV]
vThres = sqrt(Ethres*q*2/(me*(mp/me/100*mi))); %[m/s]

u = sqrt(sum(particle(4:6,:).^2,1));

figure
histogram(u)

Te3 = 0.5*sum((u > vThres).*u.^2.*particle(7,:)) / (sum(u > vThres)/me) / q

PortionAboveEnergyThreshold = sum(u > vThres) / size(particle,2);



end