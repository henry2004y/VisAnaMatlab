% Visualizing streamlines from 3D data of SWMF
%
% You can also consider streamtube function in the future. 
%
% Hongyang Zhou, hyzhou@umich.edu 08/10/2017


%% Read data One
% I found something interesting in this file.
% There is a Morton ordering from 70% to the end of the file!
% Currently I don`t need that for my Matlab visualization, but obviously
% tecplot knows how to deal with the cell indexing.
% filename = '3d__mhd_4_n00060000.dat';
% 
% data = importdata(filename,' ',23);
% 
% % the position matrix is wrong!
% xyz = data.data(:,1:3);
% rho = data.data(:,4);
% B   = data.data(:,5:7);
% hyp = data.data(:,8);
% p   = data.data(:,9);
% J   = data.data(:,10:12);

%% Read data Two
filename = '3d__mhd_4_n00060000.dat';
%filename = 'test.dat';

fid = fopen(filename,'r');

formatSpec = repmat('%f ',1,15);
N = 5632000; % Check this number in the header
data = textscan(fid,formatSpec,N,'headerlines',23);

% Still testing, not working
% for i = 1:23
%     fgetl(fid);    
% end
% data = fscanf(fid,formatSpec,[7 5632000]);

fclose(fid);

xyz = [data{1} data{2} data{3}];
%rho = data{4};
B   = [data{5} data{6} data{7}];
%hyp = data{8};
%p   = data{9};
%J   = [data{10} data{11} data{12}];

clearvars data

%% Interpolation method I
% This is too slow, and costs too much memory!
% One way to improve this is using parallel pool.
% Fx = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),B(:,1));
% Fy = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),B(:,2));
% Fz = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),B(:,3));

%xlin = linspace(plotrange(1),plotrange(2),filehead.nx(1));
%ylin = linspace(plotrange(3),plotrange(4),filehead.nx(2));
%zlin = linspace(plotrange(5),plotrange(6),filehead.nx(3));

% xlin = linspace(-3,6,30);
% ylin = linspace(-2.5,2.5,30);
% zlin = linspace(-3,3,30);
% 
% [xq,yq,zq] = meshgrid(xlin,ylin,zlin);
% bxq = Fx(xq,yq,zq);
% byq = Fy(xq,yq,zq);
% bzq = Fz(xq,yq,zq);
% 
% B   = sqrt(bxq.^2+byq.^2+bzq.^2);

%% Interpolation method II
% ranges = [-3,6,-2.5,2.5,-3,3];
% xlin = linspace(ranges(1),ranges(2),30);
% ylin = linspace(ranges(3),ranges(4),30);
% zlin = linspace(ranges(5),ranges(6),30);
% 
% [xq,yq,zq] = meshgrid(xlin,ylin,zlin);
% 
% bxq = griddata(xyz(:,1),xyz(:,2),xyz(:,3),B(:,1),xq,yq,zq);
% byq = griddata(xyz(:,1),xyz(:,2),xyz(:,3),B(:,2),xq,yq,zq);
% bzq = griddata(xyz(:,1),xyz(:,2),xyz(:,3),B(:,3),xq,yq,zq);
% 
% B   = sqrt(bxq.^2+byq.^2+bzq.^2);

%% Interpolation test: limit the region (1.6x faster than no sorting)
% tic
% ranges = [-3,6,-2.5,2.5,-3,3];
% xlin = linspace(ranges(1),ranges(2),30);
% ylin = linspace(ranges(3),ranges(4),30);
% zlin = linspace(ranges(5),ranges(6),30);
% [xq,yq,zq] = meshgrid(xlin,ylin,zlin);
% 
% % Trading speed for memory saving
% [~,ia,~] = unique(xyz * randn(size(xyz,2),1));
% xyz = xyz(ia,:);
% B = B(ia,:);
% 
% xyzIndex = xyz(:,1) > ranges(1) & xyz(:,1) < ranges(2) & ...
%            xyz(:,2) > ranges(3) & xyz(:,2) < ranges(4) & ...
%            xyz(:,3) > ranges(5) & xyz(:,3) < ranges(6) & ...
%            xyz(:,1).^2+xyz(:,2).^2+xyz(:,3).^2>1;
% 
% xyz = xyz(xyzIndex,:);
% B   = B(xyzIndex,:);
% 
% % tic
% % xyz(~xyzIndex,:) = [];
% % B(~xyzIndex,:) = [];
% % toc
% disp('sorting finished')
% 
% Fx = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),B(:,1));
% Fy = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),B(:,2));
% Fz = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),B(:,3));
% 
% disp('scatteredInterpolant created')
% 
% bxq = Fx(xq,yq,zq);
% byq = Fy(xq,yq,zq);
% bzq = Fz(xq,yq,zq);
% 
% disp('Interpolation done')
% 
% B   = sqrt(bxq.^2+byq.^2+bzq.^2);
% toc
% return

%% Parallel interpolation tests (in fact slower than pure serial code...)
tic
ranges = [-3,6,-2.5,2.5,-3,3];
% Trading speed for memory saving
[~,ia,~] = unique(xyz * randn(size(xyz,2),1));
xyz = xyz(ia,:);
B = B(ia,:);

xyzIndex = xyz(:,1) > ranges(1) & xyz(:,1) < ranges(2) & ...
           xyz(:,2) > ranges(3) & xyz(:,2) < ranges(4) & ...
           xyz(:,3) > ranges(5) & xyz(:,3) < ranges(6) & ...
           xyz(:,1).^2+xyz(:,2).^2+xyz(:,3).^2>1;

xyz = xyz(xyzIndex,:);
B = B(xyzIndex,:);

x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
Bx = B(:,1); By = B(:,2); Bz = B(:,3); 
clearvars xyz B
toc
disp('sorting finished')

%pmode reset
%pmode start local 3
%mpiprofile on
tic
parpool(3)

%spmd, profile('on','-timer','real'), end

spmd
  switch labindex
    case 1
      F = scatteredInterpolant(x,y,z,Bx);
    case 2
      F = scatteredInterpolant(x,y,z,By);
    case 3
      F = scatteredInterpolant(x,y,z,Bz);
  end
end

%spmd, p = profile('info'); profile('off'), end

% 
% Fx = parfeval(@scatteredInterpolant,1,xyz(:,1),xyz(:,2),xyz(:,3),B(:,1));
% Fy = parfeval(@scatteredInterpolant,1,xyz(:,1),xyz(:,2),xyz(:,3),B(:,2));
% Fz = parfeval(@scatteredInterpolant,1,xyz(:,1),xyz(:,2),xyz(:,3),B(:,3));
% 
% Fx = fetchOutputs(Fx);
% Fy = fetchOutputs(Fy);
% Fz = fetchOutputs(Fz);
% 
toc
% This line is extremely slow! Does this mean that most of the time is
% spent on data transfer from each worker to Matlab client?
% hyzhou: send() and poll() may save me from this (2017a)
tic
Fx = F{1}; Fy = F{2}; Fz = F{3};
toc

xlin = linspace(-3,6,30);
ylin = linspace(-2.5,2.5,30);
zlin = linspace(-3,3,30);

[xq,yq,zq] = meshgrid(xlin,ylin,zlin);
tic
bxq = Fx(xq,yq,zq);
byq = Fy(xq,yq,zq);
bzq = Fz(xq,yq,zq);
toc

B   = sqrt(bxq.^2+byq.^2+bzq.^2);

%mpiprofile('viewer',p{1})
delete(gcp)
return


%% Parallel interpolation test II (slower than serial code)
% ranges = [-3,6,-2.5,2.5,-3,3];
% xlin = linspace(ranges(1),ranges(2),30);
% ylin = linspace(ranges(3),ranges(4),30);
% zlin = linspace(ranges(5),ranges(6),30);
% [xq,yq,zq] = meshgrid(xlin,ylin,zlin);
% 
% x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
% 
% % Trading speed for memory saving
% [~,ia,~] = unique(xyz * randn(size(xyz,2),1));
% x = x(ia); y = y(ia); z = z(ia);
% B = B(ia,:);
% 
% xyzIndex = x > ranges(1) & x < ranges(2) & ...
%            y > ranges(3) & y < ranges(4) & ...
%            z > ranges(5) & z < ranges(6) & ...
%            x.^2+y.^2+z.^2>1;
% 
% x = x(xyzIndex); y = y(xyzIndex); z = z(xyzIndex);
% B = B(xyzIndex,:);
% 
% disp('sorting finished')
% tic
% F = {[],[],[]};
% parpool(3)
% parfor idx = 1:3
%    F{idx} = scatteredInterpolant(x,y,z,B(:,idx));
%    bq{idx} = F{idx}(xq,yq,zq);
% end
% toc
% disp('Interpolation done')
% 
% %B   = sqrt(bxq.^2+byq.^2+bzq.^2);
% B = sqrt(bq{1}.*bq{1}+bq{2}.*bq{2}+bq{3}.*bq{3});
% tic
% delete(gcp)
% toc
% return

%% Generating starting points
% [x,y,z] = sphere(21)
% Here is a way to generate random points on a sphere that is statistically 
% uniform. If you want a radius other than one, multiply the second line by
% the desired radius, R.
n = 100;
r = randn(3,n); % Use a large n
r = r./sqrt(sum(r.^2,1));
x = r(1,:);
y = r(2,:);
z = r(3,:);
%x = [-2.9 ;-2.8]; y = [0 ;0]; z = [1;-1];

%% Streamline with colors
% I don`t know why, but there`s a significant latency between streamline
% calculation and figure display.

tic
disp('begin plotting...')
figure(1)
p = streamlinec(stream3(xq,yq,zq,bxq,byq,bzq,x,y,z),xq, yq, zq, B);
disp('stream done')
h = colorbar; ylabel(h,'B [nT]'); hold on
disp('begin adding sphere')
sphere(21); hold off
disp('sphere added')
%colormap gray;
shading interp; light; axis equal; grid on;
xlim([-3 4.5]); ylim([-3 3]); zlim([-3 3]);
xlabel('x [$R_G$]','Interpreter','LaTex');
ylabel('y [$R_G$]','Interpreter','LaTex');
zlabel('z [$R_G$]','Interpreter','LaTex');
set(gca,'Fontsize',20);
view(3)
toc
%streamline(xq,yq,zq,bxq,byq,bzq,x,y,z)
%streamline(stream3(xq,yq,zq,bxq,byq,bzq,x,y,z))

%% Boundary (i.e. Magnetopause)
% p = findobj(p);
% pts = [0 0 0];
% for i = 1:numel(p)
%    pts = [pts ; p(i).Vertices];
% end
% k = boundary(pts);
% figure(2);
% trisurf(k,pts(:,1),pts(:,2),pts(:,3),'FaceColor','magenta'...
%    ,'Edgecolor','none','FaceAlpha',0.8,'FaceLighting','gouraud')
% %camlight
% light('Position',[-3 -1 0],'Style','infinite')