% Script for SWMF Galileo Flyby Simulation Test Comparison
% Purpose: find the best resistivity value in terms of quasi-steady state
%          magnetic field comparison for a certain flyby by selecting the 
%          smallest norm2 value.
% Prerequisite: Folder FileIO_SWMF must be included into Matlab path
% 
% Hongyang Zhou, hyzhou@umich.edu 07/31/2017

%% Parameters
flyby = 8;   % [1,2,7,8,28,29]
DoPlot = 0;  % Plot output
DoSave = 0;  % Save norm2 number
nboxfile=1;  % # of box output files
% using function griddata or scatteredInterpolant to do
% inter/extra-polation
griddataORscatteredIntepolant = 1;

%% Read observation data
flybyfile = strcat('Galileo_G',int2str(flyby),'_flyby_MAG.dat');
f = fullfile('~/Ganymede/GalileoData/galileomagdata',flybyfile);
[~,data] = read_log_data(f);

time = datetime(data(:,1:6));
xyz  = data(:,7:9);
Bobs = data(:,10:12);

%% Read/Plot simulation data for resistivity tests
norm2 = Inf(nboxfile,1);
for i=1:nboxfile
   
%filename = strcat('~/Ganymede/DivBCTtest/boxFile/box_',int2str(i),'.outs');
filename = '~/Ganymede/SeriesTest/Run78*/GM/box*.outs';
[filehead,data,filelist] = read_data(filename);

   for ipict=1:filelist.npictinfiles
   %for ipict=filelist.npictinfiles:filelist.npictinfiles
      [filehead,data] = read_data(filename,'npict',ipict);
 
      % Interpolate simulation data to observation data
      data = data.file1;
      nx = filehead.nx; nw = filehead.nw;
      
      x = data.x(:,:,:,1);
      y = data.x(:,:,:,2);
      z = data.x(:,:,:,3);
      
      w = data.w;
      % Find the correct index of variables
      bx_ = strcmpi('bx',filehead.wnames);
      by_ = strcmpi('by',filehead.wnames);
      bz_ = strcmpi('bz',filehead.wnames);
      bx = data.w(:,:,:,bx_);
      by = data.w(:,:,:,by_);
      bz = data.w(:,:,:,bz_);
      
      x  = permute(x,[2 1 3]);
      y  = permute(y,[2 1 3]);
      z  = permute(z,[2 1 3]);
      bx = permute(bx,[2 1 3]);
      by = permute(by,[2 1 3]);
      bz = permute(bz,[2 1 3]);
      
      
      if griddataORscatteredIntepolant==1
         % Using griddata will only have interpolation
         % This is about 3 times faster than scatteredInterpolant
         %    Bsim = Inf(size(xyz,1),3);
         %    Bsim(:,1) = griddata(x(:,1),x(:,2),x(:,3),w(:,5)...
         %       ,xyz(:,1),xyz(:,2),xyz(:,3),'linear');
         %    Bsim(:,2) = griddata(x(:,1),x(:,2),x(:,3),w(:,6)...
         %       ,xyz(:,1),xyz(:,2),xyz(:,3),'linear');
         %    Bsim(:,3) = griddata(x(:,1),x(:,2),x(:,3),w(:,7)...
         %       ,xyz(:,1),xyz(:,2),xyz(:,3),'linear');
         Bsim = Inf(size(xyz,1),3);
         Bsim(:,1) = interp3(x,y,z,bx,xyz(:,1),xyz(:,2),xyz(:,3));
         Bsim(:,2) = interp3(x,y,z,by,xyz(:,1),xyz(:,2),xyz(:,3));
         Bsim(:,3) = interp3(x,y,z,bz,xyz(:,1),xyz(:,2),xyz(:,3));
         
      elseif griddataORscatteredIntepolant==2
         % Using scatteredInterpolant will have both interpolation and
         % extrapolation
         Bsim = Inf(size(xyz,1),3);
         F1 = scatteredInterpolant(x(:,1),x(:,2),x(:,3),w(:,5));
         F1.Method = 'linear';
         F2 = scatteredInterpolant(x(:,1),x(:,2),x(:,3),w(:,6));
         F2.Method = 'linear';
         F3 = scatteredInterpolant(x(:,1),x(:,2),x(:,3),w(:,7));
         F3.Method = 'linear';
         
         Bsim(:,1) = F1(xyz(:,1),xyz(:,2),xyz(:,3));
         Bsim(:,2) = F2(xyz(:,1),xyz(:,2),xyz(:,3));
         Bsim(:,3) = F3(xyz(:,1),xyz(:,2),xyz(:,3));
      end
 
      % Calculate differences & save to output file
      Bdiff = Bobs - Bsim;
      % There`s a difference between 2 definitions of 2-norm.
      %norm(Bdiff(~isnan(Bdiff)))
      normNew = mean(abs(Bdiff(~isnan(Bdiff))).^2)^(1/2);
      norm2(i) = min(norm2(i),normNew);
      if norm2(i)==normNew; ipictSmall = ipict; end
      
   end
              
% Visualization
if DoPlot
   figure('Position', [100, 100, 1049, 895]);
   subplot(411); LW1 = 3; LW2 = 1.5; FS=20;
   plot(time, Bobs(:,1),time,Bsim(:,1),'k','LineWidth',LW1); %Bx
   ylabel('Bx [nT]');
   title(strcat('Galileo G',int2str(flyby),' Flyby Magnetic field'))
   set(gca,'FontSize',FS,'xMinorTick','on','LineWidth',LW2);
   subplot(412);
   plot(time, Bobs(:,2),time,Bsim(:,2),'k','LineWidth',LW1); %By
   ylabel('By [nT]');
   set(gca,'FontSize',FS,'xMinorTick','on','LineWidth',LW2);
   subplot(413);
   plot(time, Bobs(:,3),time,Bsim(:,3),'k','LineWidth',LW1); %Bz
   ylabel('Bz [nT]');
   set(gca,'FontSize',FS,'xMinorTick','on','LineWidth',LW2);
   subplot(414);
   plot(time, sqrt(Bobs(:,1).^2+Bobs(:,2).^2+Bobs(:,3).^2)...
       ,time, sqrt(Bsim(:,1).^2+Bsim(:,2).^2+Bsim(:,3).^2)...
       ,'k','LineWidth',LW1); %|B|
   ylabel('B [nT]');
   set(gca,'FontSize',FS,'xMinorTick','on','LineWidth',LW2);
end



if DoSave
   fileID = fopen('Bdiff.log','a');
   fprintf(fileID,'norm(Bdiff) %f\n',norm2);
   fclose(fileID);
end

end

% Find the index and value of the best parameter set
[val, idx] = min(norm2);

fprintf('---------------------------\n');
fprintf('Best fit output number = %d\n',idx);
fprintf('ipict, min(norm2) = %d, %f\n',ipictSmall,val);
fprintf('---------------------------\n')