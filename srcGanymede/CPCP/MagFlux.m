% Magnetic flux through the polar region.
%
% Calculate the surface integral of magnetic field of a region from a plane
% cut. 
%
% The surface integral code kernel is obtained from MathWorks Community,
% written by Richard Brown.
%
% Hongyang Zhou, hyzhou@umich.edu  01/17/2018

clear; clc
%%
flyby = 'G8'; % default is G8

Rg = 2634000; %[m], radius of Ganymede
StatusIn = 2; % Region connectivity defined by BATSRUS
StatusOut= 0; % Region connectivity defined by BATSRUS

switch flyby
   case 'G8'
      disp('G8 flyby magnetic flux calculation')
      filename= '../../newPIC/run_G8_newPIC/box_CPCP_G8_1200s.outs';
   case 'G28'
      disp('G28 flyby magnetic flux calculation')
      filename= '../../newPIC/run_G28_newPIC/box_CPCP_G28_1200s.outs';
end

[~,~,fileinfo] = read_data(filename,'verbose',false);
npict = fileinfo.npictinfiles; % # of snapshots in the file

%% Magnetic Flux through a plane surface closed region
% Right now the results seems to be wrong. I need to verify it.

Flux = Inf(npict,1);
time = Inf(npict,1);
dsample = 1/32; % resolution of the grid

% Loop over snapshots
for ipict=1:1%npict
   fprintf('ipict=%d\n',ipict);
   [filehead,data] = read_data(filename,'verbose',false,'npict',ipict);

   time(ipict) = filehead.time;
   
   %x = data.file1.x(:,:,:,1);
   %y = data.file1.x(:,:,:,2);
   status = data.file1.w(:,:,:,14);
   status(status>1) = StatusIn;
   status(status<1) = StatusOut;
   
   W = data.file1.w(:,:,:,7);
   
   [m, n] = size(W);
   
   % Average over cells
   Bz = 0.25*(W(1:m-1,1:n-1) + W(1:m-1,2:n) + W(2:m,1:n-1) + W(2:m,2:n));
   status = 0.25*(status(1:m-1,1:n-1) + status(1:m-1,2:n) + ...
      status(2:m,1:n-1) + status(2:m,2:n));
   
   % Test
   %Bz(:,:) = 1;
   
   % Eliminate the region outside the polar connected region
   Bz(status < StatusIn) = 0;

   % Total magnetic flux
   Flux(ipict) = sum(Bz(:))*dsample^2;
   
end

% Unit transformation
Flux = Flux * Rg^2 * 1e-9; %[Tm]

%% Visualization

figure;
plot(time,Flux)
xlabel('time [s]'); ylabel('Polar Flux [T\cdot m]')

%% Magnetic flux throught 3 plane cuts

filename= '../../newPIC/box_3cuts_test.out';

[filehead,data] = read_data(filename,'verbose',false,'npict',1);

data = data.file1;

dsample = 1/32; % resolution of the grid

% x = data.x(:,:,:,1);
% y = data.x(:,:,:,2);
% z = data.x(:,:,:,3);
% status = data.w(:,:,:,14);

% Loop over cuts
for icut=1:3
   status = data.w(:,:,icut,14);
   status(status>1) = StatusIn;
   status(status<1) = StatusOut;
   
   W = data.w(:,:,icut,7); % actually this is Bz
   
   [m, n] = size(W);
   
   % Average over cells
   Bz = 0.25*(W(1:m-1,1:n-1) + W(1:m-1,2:n) + W(2:m,1:n-1) + W(2:m,2:n));
   status = 0.25*(status(1:m-1,1:n-1) + status(1:m-1,2:n) + ...
      status(2:m,1:n-1) + status(2:m,2:n));
   
   % Eliminate the region outside the polar connected region
   Bz(status < StatusIn) = 0;

   % Total magnetic flux
   Flux(icut) = sum(Bz(:))*dsample^2* Rg^2 * 1e-9; %[Tm]
   
end

%%

% time derivative of magnetic flux should equal the closed line integral of
% electric field
dPhidt = gradient(Flux,1);

figure;
plot(time,dPhidt)
xlabel('time [s]'); ylabel('d\phi/dt')

%%
pi*4*Rg^2*(-160)*1e-9

Bz(Bz<0) = 1;
stats = regionprops(Bz,'area');