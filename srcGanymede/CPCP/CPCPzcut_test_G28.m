%filename = '../newPIC/G28_CPCP_test.out';
filename = '../newPIC/G28_box_z=2cut.outs';
filename = '../../newPIC/GM_G28/box_var_6_n60000_60363.outs';
%filename = '../newPIC/GM_mix2/box_var_5_n60000_60343.outs';

npict = 1;
CPCPt = Inf(npict,1);
time = Inf(npict,1);

Rg = 2634000; %[m]

for ipict=1:npict
fprintf('ipict=%d\n',ipict);
[filehead,data] = read_data(filename,'verbose',false,'npict',ipict);

time(ipict) = filehead.time;

func      = 'status';
plotmode  = 'contbar';
 
plot_data(data.file1,filehead,func,'plotmode',plotmode,'plotinterval',0.1);

%
data = data.file1;
x = data.x(:,:,:,1);
y = data.x(:,:,:,2);
z = 2*ones(size(data.x,1),size(data.x,2));
status = data.w(:,:,:,14);
status(status>1) = 2;
status(status<1) = 0;


% Get E field: E = -uxB
rho = data.w(:,:,:,1); rho = rho(:);

ux = data.w(:,:,:,2); ux = ux(:);
uy = data.w(:,:,:,3); uy = uy(:);
uz = data.w(:,:,:,4); uz = uz(:);

u = [ux uy uz];

Bx = data.w(:,:,:,5); Bx = Bx(:);
By = data.w(:,:,:,6); By = By(:);
Bz = data.w(:,:,:,7); Bz = Bz(:);
 
B = [Bx By Bz];

jx = data.w(:,:,:,11); jx = jx(:);
jy = data.w(:,:,:,12); jy = jy(:);
jz = data.w(:,:,:,13); jz = jz(:);

J = [jx jy jz];

% Calculate electron bulk velocity in Hall MHD 
% Transform the hall velocity into [km/s]
e = 1.60217662e-19; 
uex = ux - jx./rho * 1e-6 /e * 1e-3 * 1e-6;
uey = uy - jy./rho * 1e-6 /e * 1e-3 * 1e-6;
uez = uz - jz./rho * 1e-6 /e * 1e-3 * 1e-6;

ue = [uex uey uez];

%E = -cross(u,B);
E = -cross(ue,B);
E = reshape(E,size(data.w,1),size(data.w,2),3);

Ex = E(:,:,1);
Ey = E(:,:,2);
Ez = E(:,:,3);

%[curlx,curly,curlz,~] = curl(x,y,z,Ex,Ey,Ez);
% I think I understand why curl of E is not zero: E=-uxB is only valid in
% ideal MHD; in Hall MHD, I should use E=-uexB. Give it a try.
% curlx = (gradient(Ez) - gradient(Ey))*(1/32);
% curly = (gradient(Ex) - gradient(Ez))*(1/32);
% curlz = (gradient(Ey) - gradient(Ex))*(1/32);
% 
% figure;
% mesh(x,y,curlx); xlabel('x'); ylabel('y');
% figure;
% mesh(x,y,curly); xlabel('x'); ylabel('y');
% figure;
% mesh(x,y,curlz); xlabel('x'); ylabel('y');


x = x(status==2);
y = y(status==2);

% Find the boundary of magnetosphere
k = boundary(x,y,1);
figure(1); hold on
plot(x(k),y(k));

% Integrate along the boundary 
x = x(k);
y = y(k);
Ex = Ex(status==2);
Ex = Ex(k);
Ey = Ey(status==2);
Ey = Ey(k);



% Calculate cross polar cap potential (SI units!)
% Let the starting point have reference absolute potential 0
Ediff = Inf(numel(x)-1,1);
for iE=1:numel(x)-1
   Ediff(iE) = ( (Ex(iE)+Ex(iE+1))/2*(x(iE+1)-x(iE)) + ...
                 (Ey(iE)+Ey(iE+1))/2*(y(iE+1)-y(iE)) )*Rg*1e3*1e-9;
end

EPotential = cumsum([0;Ediff]);

figure;
plot(EPotential)

CPCPt(ipict) = max(EPotential) - min(EPotential);

end
% 
% figure
% plot(time,CPCPt,'*-');
% xlabel('Time [s]'); ylabel('CPCP [V]');
% set(gca,'FontSize',14,'LineWidth',1.2);

return

%%
filename = '../newPIC/G28_box_z=2cut.outs';

npict = 300;
CPCPt = Inf(npict,1);
time = Inf(npict,1);

Rg = 2634000; %[m]
dsample = 1/32;

% Background total potential drop
% G28 B= 11 78 -76 [nT], U=140 0 0 [km/s]
Uxbk = 140; Bybk = 78; Bzbk = -76;
Potential_bk = Inf(npict,1);

for ipict=1:npict
   fprintf('ipict=%d\n',ipict);
   [filehead,data] = read_data(filename,'verbose',false,'npict',ipict);
   
   time(ipict) = filehead.time;
   
   %
   x = data.file1.x(:,:,:,1);
   y = data.file1.x(:,:,:,2);
   status = data.file1.w(:,:,:,14);
   status(status>1) = 2;
   status(status<1) = 0;
   
   % Get E field: E = -uxB
   ux = data.file1.w(:,:,:,2); ux = ux(:);
   uy = data.file1.w(:,:,:,3); uy = uy(:);
   uz = data.file1.w(:,:,:,4); uz = uz(:);
   
   u = [ux uy uz];
   
   Bx = data.file1.w(:,:,:,5); Bx = Bx(:);
   By = data.file1.w(:,:,:,6); By = By(:);
   Bz = data.file1.w(:,:,:,7); Bz = Bz(:);
   
   B = [Bx By Bz];
   
   E = -cross(u,B);
   E = reshape(E,size(data.file1.w,1),size(data.file1.w,2),3);
   
   Ex = E(:,:,1);
   Ey = E(:,:,2);
   Ez = E(:,:,3);
   
   row = find(x==1,1);
   
   k = find( status(row,:)==2 );
   
   CPCPt(ipict) = -dsample * trapz( Ey(row,k) )*Rg*1e3*1e-9;
   
   Potential_bk(ipict) = -(max(y(row,k))-min(y(row,k))) * ...
      Rg * Uxbk*1e3 * (Bzbk-Bybk)*1e-9 * sind(atand(78/76));

end

figure
plot(time,CPCPt,'*-');
xlabel('Time [s]'); ylabel('CPCP [V]');
set(gca,'FontSize',14,'LineWidth',1.2);

figure
plot(time,CPCPt./Potential_bk,'*-');
xlabel('Time [s]'); ylabel('Global reconnection rate');
set(gca,'FontSize',14,'LineWidth',1.2);