% Visualizing box output from IPIC3D
%
% 1. 3D view of magnetopause
% 2. 3D view of magnetic field lines
% 3. 3D view of flux ropes
%
% I was having problem with isonormals, and it turned out to be the issue
% with meshgrid and ndgrid data structure. ndgrid structure is the
% intuitive way of thinking but this function requires meshgrid structure
% and then internally switch to ndgrid structure (such a weird thing...)
% There is a quick tranformation for the two structures.
% 
% Hongyang Zhou, hyzhou@umich.edu 09/28/2017

%{
% filename = '~/Documents/research/Ganymede/box_PIC.out';
% [filehead,data] = read_data(filename,'verbose',false);
% x = data.file1.x(:,:,:,1); 
% y = data.file1.x(:,:,:,2);
% z = data.file1.x(:,:,:,3);
% w = data.file1.w(:,:,:,7); % bz
% 
% x = linspace(-2,-1,33);
% y = linspace(-3,3,193);
% z = linspace(0.5,1,17);
% [x,y,z] = meshgrid(x,y,z);
% w = permute(w,[2 1 3]);

%[x,y,z] = meshgrid(x,y,z); % The order of x and y size is strange
%[X,Y,Z] = ndgrid(X,Y,Z);

% figure
% p = patch(isosurface(x,y,z,w,0));
% isonormals(x,y,z,w,p)
% p.FaceColor = 'red';
% p.EdgeColor = 'none';
% daspect([1 1 1])
% view(3); 
% axis tight
% camlight 
% lighting gouraud
%}

%% Read data from 3d outputs
filename = '3d*.outs'; % 3d PC outputs
[filehead,data] = read_data(filename,'verbose',false,'npict',1);

% filename = '~/Documents/research/Ganymede/box_PIC.out';
% w = data.file1.w(:,:,:,7); % bz

x = data.file1.x(:,:,:,1); 
y = data.file1.x(:,:,:,2);
z = data.file1.x(:,:,:,3);
Bx = data.file1.w(:,:,:,3);
By = data.file1.w(:,:,:,4);
Bz = data.file1.w(:,:,:,5);

% The original data is saved in ndgrid format. For streamline and
% isonormals functions, the input should be in meshgrid format.
x  = permute(x,[2 1 3]);
y  = permute(y,[2 1 3]);
z  = permute(z,[2 1 3]);
Bx = permute(Bx,[2 1 3]);
By = permute(By,[2 1 3]);
Bz = permute(Bz,[2 1 3]);

%% 3D view of magnetopause in the box region
% The magnetopause is defined by Bz==0

figure(1)
p = patch(isosurface(x,y,z,Bz,0));
%v = p.Vertices;

p.FaceColor = 'red';
p.EdgeColor = 'none';
p.FaceAlpha = 0.8;
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud
title('Dayside Magnetopause of Ganymede')
set(gca,'FontSize',16,'LineWidth',1.2)

%% 3D view of magnetic field lines
%{
% If you need refinement
% Fx = griddedInterpolant(x,y,z,Bx);
% Fy = griddedInterpolant(x,y,z,By);
% Fz = griddedInterpolant(x,y,z,Bz);
%
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
%}

%[~,v] = isosurface(x,y,z,Bz,0);
%v = v(abs(v(:,3))<0.05,:);

% Select starting pts for streamlines
[sx,sy,sz] = meshgrid(-1.72,-0.:.1:0.5,1.5);
%[sx,sy,sz] = meshgrid(-1.2,-1.6:.5:1.6,-0.8:0.125:-0.2);

% Draw field lines
hlines = streamline(x,y,z,Bx,By,Bz,sx(:),sy(:),sz(:));
%hlines = streamline(x,y,z,Bx,By,Bz,v(:,1),v(:,2),v(:,3));
%hlines = streamline(x,y,z,-Bx,-By,-Bz,v(:,1),v(:,2),v(:,3));
set(hlines,'LineWidth',2);


%% Animations
npict = 25; dpict = 1;
ipict = 1; icount = 1;
savemovie = false;

if savemovie
   figure(2)
   !mkdir -p Movie
   v = VideoWriter('Movie/movie.avi');
   v.FrameRate = 5;
   v.open
   while icount <= npict
      [filehead,data] = read_data(filename,'npict',ipict,'verbose',false);

      x  = permute(data.file1.x(:,:,:,1),[2 1 3]); 
      y  = permute(data.file1.x(:,:,:,2),[2 1 3]);
      z  = permute(data.file1.x(:,:,:,3),[2 1 3]);
      Bx = permute(data.file1.w(:,:,:,3),[2 1 3]);
      By = permute(data.file1.w(:,:,:,4),[2 1 3]);
      Bz = permute(data.file1.w(:,:,:,5),[2 1 3]);   
      
      p = patch(isosurface(x,y,z,Bz,0));
      p.FaceColor = 'red';
      p.EdgeColor = 'none';
      p.FaceAlpha = 0.8;
      daspect([1 1 1])
      view(3);
      axis([-2.5 -1 -2 2 -2 2])
      %axis square
      camlight
      lighting gouraud
      title('Dayside Magnetopause of Ganymede')
      set(gca,'FontSize',16,'LineWidth',1.2)
      
      hlines = streamline(x,y,z,Bx,By,Bz,sx(:),sy(:),sz(:));
      set(hlines,'LineWidth',2);
      
      dim = [0.13 0.13 0.1 0.046];
      str = sprintf('it=%d s',filehead.time);
      annotation('textbox',dim,'String',str,'FitBoxToText','on');
      
      frame = getframe(gcf);
      set(gca,'nextplot','replacechildren');
      clf;
      writeVideo(v,frame);     
      ipict = ipict + dpict; icount = icount + 1;
   end   
   v.close
   close(2)
end