% Local Reconnection Rate Explorer
% I need to think of a way to combine the outputs from GM and PC
%
% 1. Get magntopause position using status from GM box output for FTE;
% 2. Get Ey from PIC output;
% 3. Interpolate Ey onto the magnetopause surface.
%
% subscript 1 represents the magnetospheric side;
% subscript 2 represents the upwind side (note that there is no 
% magnetosheath at Ganymede);
%
% new idea:
% 1. find the reconnection site using E^\prime term;
% 2. create a line in the boundary normal direction through the
% reconnection site; The problem is that the strongest residual does not
% locate exactly on the fitted surface from status.
% 3. interpolate B and rho onto that line.
%
%
% Hongyang Zhou, hyzhou@umich.edu 01/23/2018
% Modified 02/22/2018

clear; clc
%% Parameters
mu0 = 4*pi*1e-7;  %[H/m]
mp  = 1.6726e-27; %[kg]
ipict = 30;
C1 = 0.95; % shrink factor for magnetosphere 
C2 = 1.1;  % expansion factor for upstream
offsetDist = 0.08;

%% Magnetopause position

filename = '~/Ganymede/newPIC/run_G8_newPIC/box_FTE_G8_1200s.outs';

[filehead,data] = read_data(filename,'verbose',false,'npict',ipict);

data = data.file1;
x = data.x(:,:,:,1); 
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);

VarIndex_ = strcmpi('status',filehead.wnames);
status = data.w(:,:,:,VarIndex_);

x3_1 = x(status==3 & x<0);
y3_1 = y(status==3 & x<0);
z3_1 = z(status==3 & x<0);

% Pick the compact boundary points coords.
s = 0.5; % compact boundary factor [0,1]
k31  = boundary(x3_1,y3_1,z3_1,s);
bc31 = unique(k31);

x3bc = x3_1(bc31);
y3bc = y3_1(bc31);
z3bc = z3_1(bc31);

rThres = 1.6; % [Rg] radius threshold for finding dayside boundary pts
xThres = -1.5;

% Find the outer boundary points
mapindex_ = (x3bc.^2 + y3bc.^2 + z3bc.^2) > rThres^2 & x3bc < xThres;
x3bc = x3bc(mapindex_);
y3bc = y3bc(mapindex_);
z3bc = z3bc(mapindex_);

% figure;
% scatter3(x3bc,y3bc,z3bc,'.');

%Fit the closed field line boundary with hypersurface

% Set up fittype and options.
ft = fittype( 'poly55' );

% Fit model to data.
[fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

% Plot fit with data.
% figure( 'Name', 'poly5' );
% h = plot( fitresult );
% legend( h, 'poly5, x=x(y,z)', 'Location', 'NorthEast' );
% % Label axes
% xlabel('x3bc [R_G]')
% ylabel('y3bc [R_G]')
% zlabel('z3bc [R_G]')
% grid on
% axis equal
% 
% xx = h.XData;
% yy = h.YData;
% zz = h.ZData;
% set(h, 'XData', zz, 'YData', xx, 'ZData', yy);
% xlim([-2 -1]); ylim([-1.1 1.1]); zlim([-0.9 0.9]);

% Generate meshgrid from fitted surface
ymin = -1.; ymax = 1.; zmin = -0.5; zmax = 0.5;
dy = 1/30; dz = dy;
[yq,zq] = meshgrid(ymin:dy:ymax,zmin:dz:zmax);
xq = fitresult(yq,zq);

%% 3d PIC outputs
filenamePC='~/Ganymede/newPIC/run_G8_newPIC/3d_fluid_35.outs';
[filehead,data] = read_data(filenamePC,'verbose',false,'npict',ipict);

data = data.file1;
x = data.x(:,:,:,1);       % [Rg]
y = data.x(:,:,:,2);
z = data.x(:,:,:,3);
Rhoe = data.w(:,:,:,1);    % [amu/cc]
Rhoi = data.w(:,:,:,2);    % [amu/cc]
Bx = data.w(:,:,:,3);      % [nT]
By = data.w(:,:,:,4);
Bz = data.w(:,:,:,5);
Ex = data.w(:,:,:,6)*1e-3; % [mV/m]
Ey = data.w(:,:,:,7)*1e-3; % [mV/m]
Ez = data.w(:,:,:,8)*1e-3; % [mV/m]
Uex = data.w(:,:,:,9);     % [km/s]
Uey = data.w(:,:,:,10);
Uez = data.w(:,:,:,11);
Uix = data.w(:,:,:,12);
Uiy = data.w(:,:,:,13);
Uiz = data.w(:,:,:,14);

% The original data is saved in ndgrid format. For streamline and
% isonormals functions, the input should be in meshgrid format.
x  = permute(x,[2 1 3]);
y  = permute(y,[2 1 3]);
z  = permute(z,[2 1 3]);
Rhoe = permute(Rhoe,[2 1 3]);
Rhoi = permute(Rhoi,[2 1 3]);
Bx = permute(Bx,[2 1 3]);
By = permute(By,[2 1 3]);
Bz = permute(Bz,[2 1 3]);
Ex = permute(Ex,[2 1 3]);
Ey = permute(Ey,[2 1 3]);
Ez = permute(Ez,[2 1 3]);
Uex = permute(Uex,[2 1 3]);
Uey = permute(Uey,[2 1 3]);
Uez = permute(Uez,[2 1 3]);
Uix = permute(Uix,[2 1 3]);
Uiy = permute(Uiy,[2 1 3]);
Uiz = permute(Uiz,[2 1 3]);


%
figure
a = slice(x,y,z,Bz,xq*C2,yq*C2,zq*C2); 
colorbar; colormap(jet);hold on
a.LineStyle = 'none';
a = slice(x,y,z,Bz,[],0,0); hold off
alpha(a,0.5); axis equal; 
a(1).LineStyle = 'none'; a(2).LineStyle = 'none';

figure
a = slice(x,y,z,Ey,xq*C2,yq*C2,zq*C2); hold on
a.LineStyle = 'none';
c = colorbar; c.Label.String = '[mV/m]'; colormap(jet); 
a = slice(x,y,z,Ey,[],0,0); hold off
alpha(a,0.5); a(1).LineStyle = 'none'; a(2).LineStyle = 'none';
axis equal
title('Dayside Magnetopause of Ganymede')
xlabel('x [R_G]'); ylabel('y [R_G]'); zlabel('z [R_G]')
set(gca,'FontSize',16,'LineWidth',1.2)

%% Boundary normal coordinate
unitz = [0 0 1]; % z-direction unit vector
% Initialize local vectors
dM = Inf(3,size(xq,1),size(xq,2)); dL = dM; dN = dM;

% Get the three local directions
for ix=1:size(xq,1)
   for iy=1:size(xq,2)
      dM(:,ix,iy) = cross([xq(ix,iy) yq(ix,iy) zq(ix,iy)],unitz);
      dL(:,ix,iy) = cross(dM(:,ix,iy),[xq(ix,iy) yq(ix,iy) zq(ix,iy)]);
      dN(:,ix,iy) = [xq(ix,iy) yq(ix,iy) zq(ix,iy)];
      
      % Normalization
      dM(:,ix,iy) = dM(:,ix,iy) / norm(dM(:,ix,iy));
      dL(:,ix,iy) = dL(:,ix,iy) / norm(dL(:,ix,iy));
      dN(:,ix,iy) = dN(:,ix,iy) / norm(dN(:,ix,iy));
   end
end

%% normal line through a point

% Electron pressure gradient + electron intertia term [mV/m]
Residual = [Ex(:) Ey(:) Ez(:)] + ...
   cross([Uex(:) Uey(:) Uez(:)],[Bx(:) By(:) Bz(:)])*1e-3;

Residual = sqrt(sum(Residual.^2,2));

Residual = reshape(Residual,size(Ex));

[M,I] = max(Residual(:));
[Ix, Iy, Iz] = ind2sub(size(Residual),I);

xResMax = x(Ix,Iy,Iz); yResMax = y(Ix,Iy,Iz); zResMax = z(Ix,Iy,Iz);

distances2 = (xq - xResMax).^2 + (yq - yResMax).^2 + (zq - zResMax).^2;

% Find the closest one
[minDistance, indexOfMin] = min(distances2(:));
[Ix, Iy] = ind2sub(size(distances2),indexOfMin);

% Find points on a line normal to the magnetopause
NormalLine = -0.1:0.01:0.2;
xLine = xq(Ix,Iy) + NormalLine*dN(1,Ix,Iy);
yLine = yq(Ix,Iy) + NormalLine*dN(2,Ix,Iy);
zLine = zq(Ix,Iy) + NormalLine*dN(3,Ix,Iy);

bxLine = interp3(x,y,z,Bx,xLine,yLine,zLine);
byLine = interp3(x,y,z,By,xLine,yLine,zLine);
bzLine = interp3(x,y,z,Bz,xLine,yLine,zLine);

RhoeLine = interp3(x,y,z,Rhoe,xLine,yLine,zLine);
RhoiLine = interp3(x,y,z,Rhoi,xLine,yLine,zLine);
RhoLine  = interp3(x,y,z,Rhoe+Rhoi,xLine,yLine,zLine);

UixLine = interp3(x,y,z,Uix,xLine,yLine,zLine);
UiyLine = interp3(x,y,z,Uiy,xLine,yLine,zLine);
UizLine = interp3(x,y,z,Uiz,xLine,yLine,zLine);
UexLine = interp3(x,y,z,Uex,xLine,yLine,zLine);
UeyLine = interp3(x,y,z,Uey,xLine,yLine,zLine);
UezLine = interp3(x,y,z,Uez,xLine,yLine,zLine);

% figure;
% plot(xLine,bxLine,xLine,byLine,xLine,bzLine,xLine,RhoLine,...
%    xLine,RhoeLine,xLine,RhoiLine);
% legend({'Bx','By','Bz','RhoTot','Rhoe','Rhoi'})

bLLine = Inf(numel(xLine),1); bMLine = bLLine; bNLine = bLLine;
UiLLine = Inf(numel(xLine),1); UiMLine = UiLLine; UiNLine = UiLLine;
UeLLine = Inf(numel(xLine),1); UeMLine = UeLLine; UeNLine = UeLLine;
% In boundary normal coordinates
% LMN: 
for iLine=1:numel(xLine)
   bLLine(iLine) = [bxLine(iLine) byLine(iLine) bzLine(iLine)]*dL(:,Ix,Iy);
   bMLine(iLine) = [bxLine(iLine) byLine(iLine) bzLine(iLine)]*dM(:,Ix,Iy);
   bNLine(iLine) = [bxLine(iLine) byLine(iLine) bzLine(iLine)]*dN(:,Ix,Iy);
   
   UiLLine(iLine) = [UixLine(iLine) UiyLine(iLine) UizLine(iLine)]*dL(:,Ix,Iy);
   UiMLine(iLine) = [UixLine(iLine) UiyLine(iLine) UizLine(iLine)]*dM(:,Ix,Iy);
   UiNLine(iLine) = [UixLine(iLine) UiyLine(iLine) UizLine(iLine)]*dN(:,Ix,Iy); 
   
   UeLLine(iLine) = [UexLine(iLine) UeyLine(iLine) UezLine(iLine)]*dL(:,Ix,Iy);
   UeMLine(iLine) = [UexLine(iLine) UeyLine(iLine) UezLine(iLine)]*dM(:,Ix,Iy);
   UeNLine(iLine) = [UexLine(iLine) UeyLine(iLine) UezLine(iLine)]*dN(:,Ix,Iy);    
end

figure;
h = plot(xLine,bLLine,'r-.',xLine,bMLine,'b-.',xLine,bNLine,'g-.'); hold on
h(1).LineWidth = 1.2; h(2).LineWidth = 1.2; h(3).LineWidth = 1.2;
h = plot(xLine,UiLLine,'r',xLine,UiMLine,'b',xLine,UiNLine,'g');
h(1).LineWidth = 1.2; h(2).LineWidth = 1.2; h(3).LineWidth = 1.2;
% plot(xLine,UeLLine,'r:',xLine,UeMLine,'b:',xLine,UeNLine,'g:');
plot(xLine,RhoLine,'k:');
xL = get(gca, 'XLim'); 
plot(xL, [0 0], '--'); hold off
xlabel('x [R_G]'); ylabel('[km/s], [nT]')
legend({'B_L','B_M','B_N','U_L','U_M','U_N','Rho_i'})
set(gca,'FontSize',16,'LineWidth',1.2)

% figure;
% plot(xLine,ULLine,xLine,UMLine,xLine,UNLine); hold on
% xL = get(gca, 'XLim'); 
% plot(xL, [0 0], '--'); hold off
% legend({'U_L','U_M','U_N'})
% set(gca,'FontSize',16,'LineWidth',1.2)

% figure;
% plot(xLine,UxLine,xLine,UyLine,xLine,UzLine);
% legend({'U_x','U_y','U_z'})
% set(gca,'FontSize',16,'LineWidth',1.2)

rho1 = RhoLine(7);
% Bx1 = bxLine(7);
% By1 = byLine(7);
% Bz1 = bzLine(7);
% 
% B1in = [Bx1 By1 Bz1]*d1(:,Ix,Iy);
% B2in = [Bx1 By1 Bz1]*d2(:,Ix,Iy);
% B1  = sqrt(Bx1.^2 + By1.^2 + Bz1.^2);

B1in = bLLine(7);
B2in = bMLine(7);
B1   = sqrt(bLLine(7).^2 + bMLine(7).^2 + bNLine(7).^2);


rho2 = RhoLine(28);
% Bx2 = bxLine(28);
% By2 = byLine(28);
% Bz2 = bzLine(28);
% B1out = [Bx2 By2 Bz2]*d1(:,Ix,Iy);
% B2out = [Bx2 By2 Bz2]*d2(:,Ix,Iy);
% B2  = sqrt(Bx2.^2 + By2.^2 + Bz2.^2);

B1out = bLLine(28);
B2out = bMLine(28);
B2   = sqrt(bLLine(28).^2 + bMLine(28).^2 + bNLine(28).^2);

theta = acosd( (B1in.*B1out + B2in.*B2out) ./ ...
   sqrt((B1in.^2+B2in.^2).*(B1out.^2+B2out.^2)) );

Epred = (1e-18/sqrt(mp*mu0)*2)*sind(theta/2)*...
   (B1*B2)^(1.5) * ( (rho1*B1 + rho2*B2)*(B1+B2) )^(-0.5); %[mV/m]


Emeas = interp3(x,y,z,Ey,xq(Ix,Iy),yq(Ix,Iy),zq(Ix,Iy)); %[mV/m]

return
%% Get the B and rho from upstreams

Rho = Rhoe + Rhoi;

% This may not be accurate enough. I can try some other approaches.
% It turns out that this is bad. 
% uix1 = interp3(x,y,z,Uix,C1*xq,C1*yq,C1*zq);
% uiy1 = interp3(x,y,z,Uiy,C1*xq,C1*yq,C1*zq);
% uiz1 = interp3(x,y,z,Uiz,C1*xq,C1*yq,C1*zq);
% 
% uix2 = interp3(x,y,z,Uix,C2*xq,C2*yq,C2*zq);
% uiy2 = interp3(x,y,z,Uiy,C2*xq,C2*yq,C2*zq);
% uiz2 = interp3(x,y,z,Uiz,C2*xq,C2*yq,C2*zq);
% 
% rho1 = interp3(x,y,z,Rho,C1*xq,C1*yq,C1*zq);
% bxv1 = interp3(x,y,z,Bx,C1*xq,C1*yq,C1*zq);
% byv1 = interp3(x,y,z,By,C1*xq,C1*yq,C1*zq);
% bzv1 = interp3(x,y,z,Bz,C1*xq,C1*yq,C1*zq);
% 
% rho2 = interp3(x,y,z,Rho,C2*xq,C2*yq,C2*zq);
% bxv2 = interp3(x,y,z,Bx,C2*xq,C2*yq,C2*zq);
% byv2 = interp3(x,y,z,By,C2*xq,C2*yq,C2*zq);
% bzv2 = interp3(x,y,z,Bz,C2*xq,C2*yq,C2*zq);

% Now use the boundary normal direction to find the upstream conditions
xq1 = xq - squeeze(offsetDist*d3(1,:,:));
yq1 = yq - squeeze(offsetDist*d3(2,:,:));
zq1 = zq - squeeze(offsetDist*d3(3,:,:));
xq2 = xq + squeeze(offsetDist*d3(1,:,:));
yq2 = yq + squeeze(offsetDist*d3(2,:,:));
zq2 = zq + squeeze(offsetDist*d3(3,:,:));

rho1 = interp3(x,y,z,Rho,xq1,yq1,zq1);
bxv1 = interp3(x,y,z,Bx,xq1,yq1,zq1);
byv1 = interp3(x,y,z,By,xq1,yq1,zq1);
bzv1 = interp3(x,y,z,Bz,xq1,yq1,zq1);

rho2 = interp3(x,y,z,Rho,xq2,yq2,zq2);
bxv2 = interp3(x,y,z,Bx,xq2,yq2,zq2);
byv2 = interp3(x,y,z,By,xq2,yq2,zq2);
bzv2 = interp3(x,y,z,Bz,xq2,yq2,zq2);

%% Cassak-Shay formula in Borovosky`s paper

% B magnitude
b1 = sqrt(bxv1.^2+byv1.^2+bzv1.^2); %[nT]
b2 = sqrt(bxv2.^2+byv2.^2+bzv2.^2); %[nT]

% Predicted electric field from hybrid Alfven velocity
% Note the units
Epred = (1e-3*(1e-9)^2*1e3/sqrt(mp))*2*(b1.*b2).^(1.5).* ...
   ( (mu0*rho1.*b1 + mu0*rho2.*b2).*(b1+b2) ).^(-0.5); %[mV/m]

% E magnitude in tangential plane
Exmeas = interp3(x,y,z,Ex,xq(:),yq(:),zq(:));
Eymeas = interp3(x,y,z,Ey,xq(:),yq(:),zq(:));
Ezmeas = interp3(x,y,z,Ez,xq(:),yq(:),zq(:));

Emeas = [Exmeas Eymeas Ezmeas]';
Etemp1 = sum(Emeas .* reshape(d1,3,[]));
Etemp2 = sum(Emeas .* reshape(d2,3,[]));
Emeas = reshape(sqrt(Etemp1.^2+Etemp2.^2),size(xq));


figure; contourf(yq,zq,Epred); colorbar
figure; contourf(yq,zq,Emeas); colorbar
figure; contourf(yq,zq,0.1*Epred-Emeas); colorbar

return

%% boundary normal coordinates for verifying Cassak-Shay formula

% Transform vectors into local coordinate system
ui2 = Inf(size(xq)); ui1 = ui2;

for iy=1:size(xq,2)
   for ix=1:size(xq,1)
      ui1(ix,iy) = [uix1(ix,iy) uiy1(ix,iy) uiz1(ix,iy)]*d3(:,ix,iy); 
      ui2(ix,iy) = [uix2(ix,iy) uiy2(ix,iy) uiz2(ix,iy)]*d3(:,ix,iy);
   end
end

% Test: this is essentially the same as the above code. However,
% surprisingly, this is slower than the former one, which is not what I
% expected.
% uTemp = [uix1(:) uiy1(:) uiz1(:)]';
% ui1 = sum(uTemp .* reshape(d3,3,[]));
% ui1 = reshape(ui1test,size(xq));

% 
Epred = b1.*b2.*(rho1.*ui1 + rho2.*ui2)./(rho1.*b1 + rho2.*b2);
Epred = Epred *1e-3; %[mV/m]

Emeas = interp3(x,y,z,Ey,xq*C2,yq*C2,zq*C2);

figure
contourf(yq,zq,Epred-Emeas,'Linestyle','none');
colorbar; colormap jet;

% figure;slice(x,y,z,Ex,-2,0,0); colorbar; axis equal
% figure;slice(x,y,z,Ey,-2,0,0); colorbar; axis equal
% figure;slice(x,y,z,Ez,-2,0,0); colorbar; axis equal


%%
sliceomatic(Ey)

%% Slices
% Based on what I learned from slices, I can use an arbitrary slice.

xslice = -2;
yslice = 0;
zslice = 0;
%[cx,cy,cz] = meshgrid(xrange,yrange,zrange);

figure;
slice(x,y,z,Ey,xslice,yslice,zslice);



%% 3D visualization



xmin = min(x(:));
xmax = max(x(:));
ymin = min(y(:));
ymax = max(y(:));
zmin = min(z(:));
zmax = max(z(:));

xrange = linspace(xmin,xmax,8);
yrange = linspace(ymin,ymax,8);
zrange = linspace(zmin,zmax,8);
[cx,cy,cz] = meshgrid(xrange,yrange,zrange);

figure
% hcone = coneplot(x,y,z,Bx,By,Bz,cx,cy,cz,10);
% 
% hcone.FaceColor = 'red';
% hcone.EdgeColor = 'none';

streamline(x,y,z,Bx,By,Bz,cx,cy,cz)
%streamtube(x,y,z,Bx,By,Bz,cx,cy,cz)

box on; grid on
view(30,40)
%daspect([2,2,1])
title('Dayside Magnetopause of Ganymede')
xlabel('x [R_G]'); ylabel('y [R_G]'); zlabel('z [R_G]')
set(gca,'FontSize',16,'LineWidth',1.2)

