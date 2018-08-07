% Script for demonstrating CPCP calculation
%
% Creating figures for demonstrating CPCP calculation from simulation.
% * CPCP calculation from BATSRUS box output;
% * Create contour plot for showing the magnetic field connectivity;
% * Create line plot for showing the potential calculation along the
% boundary curve.
%
% Hongyang Zhou, hyzhou@umich.edu 06/28/2018

clear; clc; close all
%% CPCP calculation for one snapshot
Rg = 2634000; %[m], radius of Ganymede
e  = 1.60217662e-19; % [C], electron charge

disp('G8 flyby CPCP calculation')
%filename= '~/Ganymede/newPIC/run_G8_newPIC/box_CPCP_G8_1200s.outs';
filename='~/Documents/research/Ganymede/data/box_cut_test.out';

ipict = 1;
fprintf('ipict=%d\n',ipict);
[filehead,data] = read_data(filename,'verbose',false,'npict',ipict);

data = data.file1;
time(ipict) = filehead.time;

%
x = data.x(:,:,:,1);
y = data.x(:,:,:,2);
status = data.w(:,:,:,14);
status(status>1) = 2;
status(status<1) = 0;

% Get E field: E = -uxB
rho= data.w(:,:,:,1);rho = rho(:);
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
uex = ux - jx./rho * 1e-6 /e * 1e-3 * 1e-6;
uey = uy - jy./rho * 1e-6 /e * 1e-3 * 1e-6;
uez = uz - jz./rho * 1e-6 /e * 1e-3 * 1e-6;

ue = [uex uey uez];

E = -cross(ue,B);
E = reshape(E,size(data.w,1),size(data.w,2),3);

Ex = E(:,:,1);
Ey = E(:,:,2);
Ez = E(:,:,3);


x = x(status==2);
y = y(status==2);

% Find the boundary of magnetosphere
k = boundary(x,y,1);
% figure(1); hold on
% plot(x(k),y(k));

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
EPotential = EPotential./1e3; % [kV]
[EPotentialMax,MaxIndex] = max(EPotential);
[EPotentialMin,MinIndex] = min(EPotential);

CPCPt(ipict) = EPotentialMax - EPotentialMin;

%% Contour for B topology
func      = 'status';
plotmode  = 'contf';
%
plot_data(data,filehead,func,'plotmode',plotmode,'plotinterval',0.1);

%fig1 = gcf;
axis equal
opt.XLabel = 'x [R_G]';
opt.YLabel = 'y [R_G]';
opt.FileName = 'zCutPlaneStatus.eps';

xlabel(opt.XLabel); ylabel(opt.YLabel);
set(gca,'LineWidth',1.2)

%% Line plot for potential calculation
opt.XLabel = 'Position along the curve'; % xlabel
opt.YLabel = 'Electric Potential [kV]'; %ylabel
opt.Colors = [ % two colors for two data set
    0,   0.5,  0;
    1,   0,    0;
    ];
%opt.XTick = '';
%opt.Legend.Visible = 'off';
%opt.BoxDim = [7, 3]; %[width, height]
% Save? comment the following line if you do not want to save
opt.FileName = 'PotentialIntegral.eps'; 

fig2 = figure;
plot(1:MinIndex,EPotential(1:MinIndex),MinIndex:numel(EPotential),...
   EPotential(MinIndex:end));

% Create textarrow
annotation(fig2,'textarrow',[0.462365591397849 0.462365591397849],...
   [0.315326530612245 0.251700680272109],'String',{'(1.3,2.7)'});

% Create textarrow
annotation(fig2,'textarrow',[0.78494623655914 0.779569892473118],...
   [0.662265306122449 0.73469387755102],'String',{'(1,-2.5)'});

% Create textarrow
annotation(fig2,'textarrow',[0.173835125448029 0.134408602150538],...
   [0.789115646258503 0.755102040816326],'String',{'(1,-2.5)'});

% Create line
annotation(fig2,'line',[0.453405017921147 0.53405017921147],...
   [0.237 0.237],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',1.2);

% Create line
annotation(fig2,'line',[0.734767025089606 0.806451612903226],...
   [0.755 0.755],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',1.2);

ax = gca;
ax.XTickLabel = '';

% create the plot
setPlotProp(opt);

ax.Legend.Visible = 'off';


