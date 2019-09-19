% Output data processing for BATSRUS
%
% I cannot draw B/j. Add these functions!
%
% Hongyang Zhou, hyzhou@umich.edu 06/14/2017
%--------------------------------------------------------------------------


%% Test of log plot
filename = 'logTest.log';
% I found that the filetype argument is useless here! Remove it!
[filehead,logdata] = read_data(filename,'filetype','log');

[filehead,logdata] = read_data(filename);
plot_log_data(logdata,filehead,'pmin','plotmode','line');

%% Test of 2D plot, generalized coordinates
filename = '~/SWMF/SWMF/GM/BATSRUS/run_test/GM/IO2/y*';

[filehead,data] = read_data(filename);

plot_data(data.file1,filehead,'p','plotmode','contbar')

npict = 61;        % index of snapshot to be read  
filename = '../FileIO_SWMF/y=0_var_1_n0_60000.outs ../FileIO_SWMF/z=0_var_2_n0_60000.outs';
[filehead,data] = read_data(filename,'npict',npict);

plot_data(data.file1,filehead(1),'p','plotmode','grid');
plot_data(data.file2,filehead(2),'p','plotmode','grid');

func = 'p';
plotmode='contbar';
plot_data(data.file1,filehead(1),func,'plotmode',plotmode,...
   'plotinterval',0.2);

plotrange = [-4 4 -4 4];
plot_data(data.file1,filehead(1),func...
   ,'plotmode','trimesh','plotrange',plotrange)
 
plot_data(data.file1,filehead.file1,func...
   ,'plotmode','trisurf','plotrange',plotrange)

plot_data(data.file1,filehead.file1,func...
   ,'plotmode','meshbar','plotrange',plotrange,'plotinterval',0.05)


filename='y=0*.out';
[filehead,data] = read_data(filename);
 
func = 'jy ux;uz';
plotmode = 'contbar streamover';
plotrange = [-8 8 -8 8];
 
plot_data(data.file1,filehead(1),func...
   ,'plotmode',plotmode,'plotrange',plotrange,'plotinterval',0.05)
 
rectangle('Position',[-.5,-.5,1,1],'Curvature',[1,1]...
   ,'FaceColor',[.6 .6 .6])
viscircles([0 0],1,'color',[.4 .2 1]);

 
filename='z=0_var_2_n00060000.out';
filename='y=0_var_1_n00060000.out';
[filehead,data] = read_data(filename,'npict',1,'verbose',false);
 
func='jx ux;uy';
plotmode='contbar streamover';
plotrange = [-4 4 -4 4];
plot_data(data.file1,filehead(1),func,'plotmode',plotmode,...
   'plotrange',plotrange,'plotinterval',0.05)
hold on

rectangle('Position',[-2.8 -3 (-1.125+2.8) 6],'EdgeColor','r',...
  'LineWidth',2)

rectangle('Position',[-.5,-.5,1,1],'Curvature',[1,1]...
   ,'FaceColor',[.6 .6 .6]);
viscircles([0 0],1,'color',[.4 .2 1]);
 
hold off

func = 'p bx;by';
plotmode = 'contbar streamover';
plotrange = [-4 4 -4 4];

plot_data(data.file1,filehead(1),func...
  ,'plotmode',plotmode,'plotrange',plotrange,'plotinterval',0.05)

rectangle('Position',[-.5,-.5,1,1],'Curvature',[1,1]...
   ,'FaceColor',[.6 .6 .6]);
viscircles([0 0],1,'color',[.4 .2 1]);

func = 'jr jx;jy';
plotmode = 'contbar quiverover';
plotrange = [-4 4 -4 4];

func='jx;jy';
plotmode= 'quiver';

plot_data(data.file1,filehead(1),func,...
   'plotmode',plotmode,'plotrange',plotrange,'plotinterval',0.05)
viscircles([0 0],1,'color',[.6 .4 .8],'LineWidth',1);
xlabel('x'); ylabel('y'); 
%ylabel('z');
axis equal

%% Jr contour plot
filename='../FileIO_SWMF/z=0_var_2_n00060000.out';
%filename='../FileIO_SWMF/y=0_var_1_n00060000.out';
[filehead,data] = read_data(filename);

x = data.file1.x; w = data.file1.w;
j = sqrt(w(:,:,10).^2+w(:,:,11).^2);
jr = j.* x(:,:,2) ./ sqrt(x(:,:,1).^2+x(:,:,2).^2);

plotrange = [-1.05 1.05 -1.05 1.05]; plotinterval=0.01;
X = reshape(x(:,:,1),[],1);
Y = reshape(x(:,:,2),[],1);

xyIndex = X > plotrange(1) & X < plotrange(2) & ...
   Y > plotrange(3) & Y < plotrange(4);
X = X(xyIndex);
Y = Y(xyIndex);
jr = jr(xyIndex);

F = scatteredInterpolant(X,Y,jr);
[xq,yq] = meshgrid(plotrange(1):plotinterval:plotrange(2)...
   ,plotrange(3):plotinterval:plotrange(4));
vq = F(xq,yq);

figure(1)
contourf(xq,yq,vq,'Edgecolor','none'); c = colorbar;

rectangle('Position',[-.5,-.5,1,1],'Curvature',[1,1]...
   ,'FaceColor',[.6 .6 .6]);
viscircles([0 0],1,'color',[.6 .4 .8],'LineWidth',1);

xlabel('x'); %ylabel('y'); 
ylabel('z');
ylabel(c,'$\mu A/m^2$','Interpreter','LateX')
title('$j_r$ in meridional cut','Interpreter','LateX')
set(gca,'Fontsize',16); axis equal

figure(2)
hist(jr,1000);
xlabel('$j_r$ [$\mu A/m^2$]','Interpreter','LateX')
ylabel('Count')
set(gca,'Fontsize',16);

%% Test of animation
func = 'P';
plotmode = 'contbar';
plotrange = [-4 4 -4 4];

animate_data(filename,func,'plotmode',plotmode,'plotrange',plotrange,...
   'firstpict',1,'dpict',1,'npictmax',61,...
   'plotinterval',0.05,'savemovie',true)

%animate_data(data,filehead,func,'savemovie','png')

%% Test of reading 3D box output, Cartesian coordinates
filename = 'box*.outs';
[filehead,data] = read_data(filename);

%% Test of reading flyby trajectory data & compare with simulation outputs
[filehead,data] = read_log_data('Galileo_G8_flyby_MAG.dat');

time = datetime(data(:,1:6));
xyz  = data(:,7:9);
B_obs= data(:,10:12);
%plot(time,B_obs(:,3)); %Bz

filename = 'box*.outs';
[filehead,data] = read_data(filename);

%% Test of 2D simulation output, Cartesian coordinates
filename = '../../../research/Ganymede/newPIC/GM_G8/box_var_6_n60000_60321.outs';
npict = 1;

[filehead,data] = read_data(filename,'npict',npict);

plot_data(data.file1,filehead(1),'p','plotmode','contbar')

%% Test of coarse grid faceBC, Ganymede
filename = 'y=0_var_1_n00000000.out';
[filehead,data] = read_data(filename);

plot_data(data.file1,filehead.file1,'p','plotmode','grid');

plotrange = [-4 4 -4 4];
plot_data(data.file1,filehead.file1,'ux','plotmode','contbar',...
   'plotrange',plotrange,'plotinterval',0.05)

%% Test of tristream (slow!!!)
filename='../FileIO_SWMF/y=0_var_1_n00060000.out';
[filehead,data] = read_data(filename);

%func = 'jy ux;uz';
%plotmode = 'contbar streamover';
func = 'jy';
plotmode ='tristream';
plotrange = [-4 4 -4 4];

plot_data(data.file1,filehead(1),func...
  ,'plotmode',plotmode,'plotrange',plotrange)

rectangle('Position',[-.5,-.5,1,1],'Curvature',[1,1]...
   ,'FaceColor',[.6 .6 .6])
viscircles([0 0],1,'color',[.4 .2 1]);

%% Test of reading IPIC3D output
filename = 'y=0_fluid_region0_1_n0_35781.outs';

[filehead,data] = read_data(filename,'npict',60);

plot_data(data.file1,filehead,'ps0','plotmode','contbar')


%plotrange = [-4 4 -4 4];

animate_data(filename,'ps0',...
   'firstpict',1,'dpict',1,'npictmax',31,...
   'plotinterval',0.05,'savemovie',true)


%info = h5info(filename)
%h5disp(filename)


%%
% I think a better way next is to make a GUI for plotting.
% Play and have fun!
