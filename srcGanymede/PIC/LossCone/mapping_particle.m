function mappingParticle(particle,angle,phi1,theta1,Bsurf,B_P)
%MAPPINGPARTICLE Test of plotting on a sphere without Mapping ToolBox.
%
%INPUTS:
% particle: positions, velocities and weights, [7,nP]
% angle: pitch angle
% phi1:
% theta1:
% Bsurf:
% B_P
%

% Mapping particles to the surface
nP = size(particle,2);
vPar = Inf(nP,1); vPerp = Inf(nP,1);
% Calculate the perpendicular and parallel velocity
for iP=1:nP
   vTot = sqrt(particle(4,iP)^2 +  particle(5,iP)^2 + particle(6,iP)^2);
   vPar(iP) = vTot * cosd(angle(iP));
   vPerp(iP) = vTot * sind(angle(iP));
end

vPerp = sqrt(Bsurf./B_P) .* vPerp;
vPar  = sqrt(sum(particle(4:6,:).^2,1)' - vPerp.^2);

% Create surface mesh
%theta = -90:5:90;
%phi   = 5:5:360;
%[theta,phi] = meshgrid(theta,phi);

% Original locations of particles inside the loss cone
figure
scatter3(particle(1,:),particle(2,:),particle(3,:),'.');
axis equal

% Distribution on the surface
figure
hist3([phi1,theta1],'CdataMode','auto','Nbins',[30 30]);
%xlim([0,360]); ylim([-90,90])
xlabel('longitude'); ylabel('latitude')

% Calculate average velocities by bins
%[N,edges] = histcounts()
bins = 30;
phiMin = min(phi1); phiMax = max(phi1);
thetaMin = min(theta1); thetaMax = max(theta1);
dPhi = (phiMax - phiMin)/(bins - 1); 
dTheta = (thetaMax - thetaMin)/(bins - 1);
subs = [round((phi1-phiMin)/dPhi)+1 round((theta1-thetaMin)/dTheta)+1]; 

% vParAll = accumarray(subs,vPar.*particle(7,:)') ./ ...
%    accumarray(subs,particle(7,:)');


% Generate a sphere with lines
R = 1;
latspacing = 10; 
lonspacing = 20; 
% lines of longitude: 
[lon1,lat1] = meshgrid(-180:lonspacing:180,linspace(-90,90,300)); 
[x1,y1,z1] = sph2cart(lon1*pi/180,lat1*pi/180,R); 
figure
plot3(x1,y1,z1,'-','color',0.5*[1 1 1])
hold on
% lines of latitude: 
[lat2,lon2] = meshgrid(-90:latspacing:90,linspace(-180,180,300)); 
[x2,y2,z2] = sph2cart(lon2*pi/180,lat2*pi/180,R); 
plot3(x2,y2,z2,'-','color',0.5*[1 1 1])
axis equal tight off

[X,Y,Z] = sphere(100); 
surf(X*R*.99,Y*R*.99,Z*R*.99,'facecolor','w','edgecolor','none')

% contour mapping
theta = thetaMin:dTheta:thetaMax;
phi = phiMin:dPhi:phiMax;
[Phi,Theta] = meshgrid(phi,theta);

[X,Y,Z] = sph2cart(Phi./180*pi,Theta./180*pi,1);
[n,c] = hist3([phi1,theta1],'CdataMode','auto','Nbins',[30 30]);
%surf(Phi,Theta,n');
surf(X,Y,Z,n','EdgeColor','none')
% res = 21;
% lambda = linspace(-pi,pi,res);
% theta = linspace(-pi/2,pi/2,ceil(res/2));
% [L,T] = meshgrid(lambda,theta);
% [X,Y,Z] = sph2cart(L,T,1);
% surf(X,Y,Z,'FaceColor','none');
% axis equal

% figure;
% contourf(thetaMin:dTheta:thetaMax,phiMin:dPhi:phiMax,vParAll); colorbar

end

