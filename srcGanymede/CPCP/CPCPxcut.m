% Calculating cross polar cap potential for G8 flyby from x=0 cuts.
% Procedure:
% 1. Calculate E = -uxB from MHD output. It should be more accurate to use 
%    ue x B, where ue = u - nJ/e. What I found from the simulation is that 
%    this term is relatively small, so it is like a correction term.
% 2. Use status to find the region of half open field lines that have one
%    end connected to the polar cap; specifically in x=0 cut, I need a line
%    at a certain height within the region.
% 3. Get electric potential by integrating the electric field along the 
%    line.
%
% This script works for G8 and G28 Galileo flybys. 
%
% The interpolation is slow because it is using scattered interpolation.
% For x=0 cut in spherical coordinates, the data are arranged like
% unstructured mesh (DxSavePlot=-1).
%
% As we discussed heavily, using a line integral in the x=0 cut is not the
% proper way of calculating the electric field integral. A better way would
% be integrating along the magnetopause edges in a plane cut above/below 
% the moon. This is done in *CPCP_zcut.m*.
%
% Hongyang Zhou, hyzhou@umich.edu 11/02/2017
%
% Modified 01/04/2018
% 1. One script for both G8 and G28 flybys.
% 2. Reorganize each part of the code.
%
% Modified 06/29/2018
% Added comments.

clear; clc
%% Parameters
% All the following sections need parameters here.
flyby = 'G28'; % default is G8

Rg = 2634000; %[m], radius of Ganymede

switch flyby
   case 'G8'
      disp('G8 flyby CPCP calculation')
      filename = '~/Ganymede/newPIC/run_G8_newPIC/x=0_var_1_n60000_247407.outs'; % 2d GM outputs
      plotrange = [-5 5 -8 8];
   case 'G28'
      disp('G28 flyby CPCP calculation')
      filename = '~/Ganymede/newPIC/run_G28_newPIC/x=0_var_1_n60000_263035.outs'; % 2d GM outputs
      plotrange = [-8 8 -8 8];
end

%% Plots for visual reference
[filehead,data] = read_data(filename,'verbose',false);

data = data.file1;

ux = data.w(:,:,2);
uy = data.w(:,:,3);
uz = data.w(:,:,4);
bx = data.w(:,:,5);
by = data.w(:,:,6);
bz = data.w(:,:,7);

func      = 'status';
plotmode  = 'contbar';

plot_data(data,filehead,func,'plotrange',plotrange,...
   'plotmode',plotmode,'plotinterval',0.1);
axis equal; %colormap(jet)
set(gca,'Xdir','reverse')

%plot_data(data,filehead,'status by;bz','plotrange',[-10 10 -10 10],...
%   'plotmode','contbar streamover','plotinterval',0.1);

% func = 'ux by;bz';
% plotmode = 'contbar streamover';
% plot_data(data,filehead,func,'plotrange',plotrange,...
%    'plotmode',plotmode,'plotinterval',0.05);
% axis equal

%% Test line integral of electric field
% The line selected intercept z axis at 3[Rg] and has the slope of 11/77.
% This gives a rough estimation of the "potential".

status = data.file1.w(:,:,14);
y = data.file1.x(:,:,1); 
z = data.file1.x(:,:,2);

status(status==-1.5) = -1;
status(status==0.5)  = 1;
status(status==2.5)  = 2;

% Create a line
dsample = 0.01;

switch flyby
   case 'G8'
      % I am not sure about the sign of a
      a = 1/7; b=3;
      yq = -3:dsample:3;
   case 'G28'
      % I am not sure about the sign of a
      a = 78/76; b=2;
      yq = -5:dsample:2;      
end

zq = a*yq + b;

hold on; plot(yq,zq,'k'); hold off

uxq = griddata(y,z,ux,yq,zq);
uyq = griddata(y,z,uy,yq,zq);
uzq = griddata(y,z,uz,yq,zq);
bxq = griddata(y,z,bx,yq,zq);
byq = griddata(y,z,by,yq,zq);
bzq = griddata(y,z,bz,yq,zq);
statusq = griddata(y,z,status,yq,zq,'nearest');

% find the index of points inside the magnetosphere
index = find( statusq == 2 );

% Background total potential drop
switch flyby
   case 'G8'
      % G8 B= -10 -6 -86 [nT], U=140 0 0 [km/s], Width=5 Rg
      Width = 35*sind(8.13); Uxbk = 140; Bybk = -6; Bzbk = -86;
   case 'G28'
      % G28 B= -6 78 -76 [nT], U=140 0 0 [km/s], Width=5 Rg
      Width = 5.0297; Uxbk = 140; Bybk = 78; Bzbk = -76;

end

Potential_bk = -Width * Rg * Uxbk*1e3 * (Bzbk - Bybk)*1e-9

% Integrate along the line (initially this is negative, and I add a - sign
% at the beginning for plotting and other calculations)
CPCP = -dsample * Rg*1e3*1e-9 * ...
      (trapz( uxq(index).*bzq(index)-uzq(index).*bxq(index) + ...
       -uxq(index).*byq(index)+uyq(index).*bxq(index)))
    
fprintf('Global reconnection ratio=%f\n',CPCP/Potential_bk)
    
%% Check CPCP from one snapshot as a function of z
% Ideally, the value should be the same no matter what z it is. However, in
% practice I do see a difference.

b = linspace(2,5,10);
CPCP = zeros(10,1);

for ib=1:numel(b)
   zq = a*yq + b(ib);
   
   uxq = griddata(y,z,ux,yq,zq);
   uyq = griddata(y,z,uy,yq,zq);
   uzq = griddata(y,z,uz,yq,zq);
   bxq = griddata(y,z,bx,yq,zq);
   byq = griddata(y,z,by,yq,zq);
   bzq = griddata(y,z,bz,yq,zq);
   statusq = griddata(y,z,status,yq,zq,'nearest');

   % find the index of points inside the magnetosphere
   index = find( statusq == 2 );
   
   CPCP(ib) = -dsample * Rg*1e3*1e-9 * ...
      (trapz( uxq(index).*bzq(index)-uzq(index).*bxq(index)) + ...
       trapz( -uxq(index).*byq(index)+uyq(index).*bxq(index)));
end
   
figure(2);
plot(b,CPCP,'*-');
xlabel('b [R_G]'); ylabel('CPCP [V]');
set(gca,'FontSize',14,'LineWidth',1.2);

%% Tranform into regular Cartesian grid
% Interpolate onto regular grids and integrate along horizontal lines.

% The average of 0 and 2 is 1, so the boundary of 1 is not clear due to
% interpolation of the cut between cells. If it is a truly 3D output, then
% there`s no such problem.
status = data.file1.w(:,:,14);
y = data.file1.x(:,:,1); 
z = data.file1.x(:,:,2);

% It seems that status==0.5 is the right on the magnetopause. below x=0
% status=2.5 is the boundary between closed field line and half open ones.
%histogram(status)

status(status==-1.5) = -1;
status(status==0.5)  = 1;
status(status==2.5)  = 2;

dsample = 0.01;
[ymesh,zmesh] = meshgrid(-6:dsample:6,2:dsample:6);

uxq = griddata(y,z,ux,ymesh,zmesh);
uyq = griddata(y,z,uy,ymesh,zmesh);
uzq = griddata(y,z,uz,ymesh,zmesh);
bxq = griddata(y,z,bx,ymesh,zmesh);
byq = griddata(y,z,by,ymesh,zmesh);
bzq = griddata(y,z,bz,ymesh,zmesh);
statusq = griddata(y,z,status,ymesh,zmesh,'nearest');

CPCP = zeros(11,1);
zline = linspace(2,5,11);

for iz=1:numel(zline)
   % find the row index for the horizontal line
   row = find(zmesh==zline(iz),1); 
   
   k = find( statusq(row,:) == 2 );
   
   CPCP(iz) = dsample * ...
      trapz( -uxq(row,k).*bzq(row,k)+uzq(row,k).*bxq(row,k))*Rg*1e3*1e-9;
end

figure(3);
plot(zline,CPCP,'*-');
xlabel('z [R_G]'); ylabel('CPCP [V]');
set(gca,'FontSize',14,'LineWidth',1.2);

%% Plot CPCP calculation with b line integral as a function of time
npict = 600;
CPCPt = zeros(npict,1);
Rg = 2634000; %[m]

% Create a line
dsample = 0.01;

switch flyby
   case 'G8'
      filename='../../newPIC/old/x=0_G8_600s.outs';
      % I am not sure about the sign of a
      a = 1/7; b=3;
      yq = -3:dsample:3;
   case 'G28'
      filename = '../../newPIC/old/x=0_G28_600s.outs'; % 2d GM outputs
      % I am not sure about the sign of a
      a = 78/76; b=2;
      yq = -5:dsample:2;      
end

zq = a*yq + b;

for ipict=1:npict
   [~,data] = read_data(filename,'verbose',false,'npict',ipict);
   fprintf('ipict=%d',ipict);

   status = data.file1.w(:,:,14);
   y = data.file1.x(:,:,1);
   z = data.file1.x(:,:,2);
   ux = data.file1.w(:,:,2);
   uy = data.file1.w(:,:,3);
   uz = data.file1.w(:,:,4);
   bx = data.file1.w(:,:,5);
   by = data.file1.w(:,:,6);
   bz = data.file1.w(:,:,7);
    
   status(status==-1.5) = -1;
   status(status==0.5)  = 1;
   status(status==2.5)  = 2;
 
   uxq = griddata(y,z,ux,yq,zq);
   uyq = griddata(y,z,uy,yq,zq);
   uzq = griddata(y,z,uz,yq,zq);
   bxq = griddata(y,z,bx,yq,zq);
   byq = griddata(y,z,by,yq,zq);
   bzq = griddata(y,z,bz,yq,zq);
   statusq = griddata(y,z,status,yq,zq,'nearest');
   
   % find the index of points inside the magnetosphere
   index = find( statusq == 2 );   
   
   % Integrate along the line (initially this is negative, and I add a - sign
   % at the beginning for plotting and other calculations)
   CPCPt(ipict) = -dsample * Rg*1e3*1e-9 * ...
         (trapz( uxq(index).*bzq(index)-uzq(index).*bxq(index) + ...
         -uxq(index).*byq(index)+uyq(index).*bxq(index)));   
   
end

figure(4);
plot(1:npict,CPCPt,'*-');
xlabel('Time [s]'); ylabel('CPCP [V]');
set(gca,'FontSize',14,'LineWidth',1.2);

%% FFT
% Performing fft on CPCP to get period information
CPCPfft = fft(CPCPt);
CPCPfft(1) = [];
n = length(CPCPfft);
% power of first half of transform data
power = abs(CPCPfft(1:floor(n/2))).^2; 
maxfreq = 1/2;                   % maximum frequency
freq = (1:n/2)/(n/2)*maxfreq;    % equally spaced frequency grid
figure
plot(freq,power)
xlabel('Cycles/Second')
ylabel('Power')

period = 1./freq;
figure;
plot(period,power);
xlim([0 50]); % zoom in on max power
xlabel('Seconds/Cycle')
ylabel('Power')

%% Plot CPCP calculation with z=2 line integral as a function of time

npict = 600;
CPCPt = zeros(npict,1);
time = zeros(npict,1);
Rg = 2634000; %[m]

dsample = 0.02;
[ymesh,zmesh] = meshgrid(-4:dsample:4,2:dsample:2);

switch flyby
   case 'G8'
      filename = '../../newPIC/old/x=0_G8_600s.outs';
   case 'G28'
      filename = '../../newPIC/old/x=0_G28_600s.outs';  
end

for ipict=1:npict
   [filehead,data] = read_data(filename,'verbose',false,'npict',ipict);
   time(ipict) = filehead.time;
   fprintf('ipict=%d\n',ipict);
   
   status = data.file1.w(:,:,14);
   ux = data.file1.w(:,:,2);
   uz = data.file1.w(:,:,4);
   bx = data.file1.w(:,:,5);
   bz = data.file1.w(:,:,7);
   y = data.file1.x(:,:,1);
   z = data.file1.x(:,:,2);
   
   status(status==-1.5) = -1;
   status(status==0.5)  = 1;
   status(status==2.5)  = 2;
   
   uxq = griddata(y,z,ux,ymesh,zmesh);
   uzq = griddata(y,z,uz,ymesh,zmesh);
   byq = griddata(y,z,bx,ymesh,zmesh);
   bzq = griddata(y,z,bz,ymesh,zmesh);
   statusq = griddata(y,z,status,ymesh,zmesh,'nearest');
   
   row = find(zmesh==2,1);
   
   k = find( statusq(row,:) == 2 );
   
   CPCPt(ipict) = dsample * ...
      trapz( -uxq(row,k).*bzq(row,k)+uzq(row,k).*byq(row,k))*Rg*1e3*1e-9;
end

figure(5);
plot(time,CPCPt,'*-');
xlabel('Time [s]'); ylabel('CPCP [V]');
set(gca,'FontSize',14,'LineWidth',1.2);

%% Old method (not used!)
%{
%% Find the positions of the status boundary
% The average of 0 and 2 is 1, so the boundary of 1 is not clear due to
% intepolation of the cut. If it is a truly 3D output, then there`s no such
% problem.
% status = data.file1.w(:,:,14);
% x = data.file1.x(:,:,1); 
% z = data.file1.x(:,:,2);
% 
% % It seems that status==0.5 is the right on the magnetopause. below x=0
% % status=2.5 is the boundary between closed field line and half open ones.
% %histogram(status)
% 
% status(status==-1.5) = -1;
% status(status==0.5)  = 1;
% status(status==2.5)  = 2;
% 
% data.file1.w(:,:,14) = status;
% 
% % Pick the pts on the boundary
% x1 = x(status==2 & z>0);
% z1 = z(status==2 & z>0);
% k1 = boundary(x1,z1,0.5); % The factor can be adjusted between [0,1].
% 
% plot_data(data.file1,filehead,func,'plotrange',plotrange,...
%    'plotinterval',0.1);
% hold on
% plot(x1(k1),z1(k1),'k','linewidth',1.2); 
% hold off; axis equal
% 
% % Pick the pts on the boundary
% x2 = x(status==1 & z<0);
% z2 = z(status==1 & z<0);
% k2 = boundary(x2,z2,0.5); % The factor can be adjusted between [0,1].
% 
% plot_data(data.file1,filehead,func,'plotrange',plotrange,...
%    'plotinterval',0.1);
% hold on
% plot(x2(k2),z2(k2),'k','linewidth',1.2); 
% hold off; axis equal

%% Calculate CPCP 
% % Integrate along a horizontal line along y direction
% % y = 2
% xy1 = [-2.33 1.936]; xy2 = [2.366 2.067];
% [xmesh,ymesh] = meshgrid(-2.33:0.001:2.366,1.936:0.001:2.067);
% % y = 2.5
% %xy1 = [-2.2285 2.491]; xy2 = [2.359 2.508];
% %[xmesh,ymesh] = meshgrid(-2.229:0.001:2.359,2.491:0.001:2.508);
% % y = 3
% %xy1 = [-2.153 3.097]; xy2 = [2.308 3.072];
% %[xmesh,ymesh] = meshgrid(-2.153:0.001:2.308,3.072:0.001:3.097);
% % y = 4
% %xy1 = [-1.827 4.055]; xy2 = [2.111 4.008];
% %[xmesh,ymesh] = meshgrid(-1.827:0.001:2.111,4.008:0.001:4.055);
% % y = 5
% %xy1 = [-1.421 4.899]; xy2 = [1.721 5.003];
% %[xmesh,ymesh] = meshgrid(-1.421:0.001:1.721,4.899:0.001:5.003);
% % y = 6
% %xy1 = [-0.6751 6.087]; xy2 = [1.121 6.021];
% %[xmesh,ymesh] = meshgrid(-0.676:0.001:1.121,6.021:0.001:6.087);
% 
% % y = -2
% %xy1 = [-2.604 -2.058]; xy2 = [2.489 -1.918];
% %[xmesh,ymesh] = meshgrid(-2.604:0.001:2.489,-2.058:0.001:-1.918);
% % y = -3
% %xy1 = [-2.707 -3.024]; xy2 = [2.492 -2.925];
% %[xmesh,ymesh] = meshgrid(-2.707:0.001:2.492,-3.024:0.001:-2.925);
% 
% 
% t = linspace(0,1,250)';
% xyseg = (1-t)*xy1 + t*xy2;
% 
% 
% 
% uxq = griddata(x,z,ux,xmesh,ymesh);
% uzq = griddata(x,z,uz,xmesh,ymesh);
% bxq = griddata(x,z,bx,xmesh,ymesh);
% bzq = griddata(x,z,bz,xmesh,ymesh);
% 
% uxseg = interp2(xmesh,ymesh,uxq,xyseg(:,1),xyseg(:,2));
% uzseg = interp2(xmesh,ymesh,uzq,xyseg(:,1),xyseg(:,2));
% bxseg = interp2(xmesh,ymesh,bxq,xyseg(:,1),xyseg(:,2));
% bzseg = interp2(xmesh,ymesh,bzq,xyseg(:,1),xyseg(:,2));
% 
% d = cumsum([0;sqrt(sum(diff(xyseg).^2,2))]);
% % Calculate cross polar cap potential (SI units!)
% CPCP = trapz(d,-uxseg.*bzseg+uzseg.*bxseg)*Rg*1e3*1e-9

%}


