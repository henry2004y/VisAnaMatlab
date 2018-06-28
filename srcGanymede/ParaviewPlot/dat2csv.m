% Create traj.csv from .dat format
%
% Hongyang Zhou, hyzhou@umich.edu 03/24/2018

%%
flyby = 28;
flybyfile = strcat('Galileo_G',int2str(flyby),'_flyby_MAG.dat');
f = fullfile('~/Documents/research/Ganymede/Galileo',flybyfile);
[~,data] = read_log_data(f);

xyz  = data(:,7:9);

xyz( mod(1:size(xyz,1),25)>0,:) = [];

filename = 'test_G28.csv';
fid = fopen(filename, 'wt');
fprintf(fid, 'x,y,z\n');
dlmwrite(filename,xyz,'-append')