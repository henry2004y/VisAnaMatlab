%% An ugly try of grid cut plot for Ganymede simulation

npict = 61;        % index of snapshot to be read
filename = 'y=0_var_1_n0_60000.outs z=0_var_2_n0_60000.outs';
[filehead,data] = read_data(filename,'npict',npict);

nxz = filehead.file1.nx; nxy = filehead.file2.nx;
y = zeros(nxz(1),nxz(2)); z = zeros(nxy(1),nxy(2));
figure
scatter3(data.file1.x(:,1,1),y,data.file1.x(:,1,2),'.','LineWidth',1);
hold on
scatter3(data.file2.x(:,1,1),data.file2.x(:,1,2),z,'.','LineWidth',1);
%hold off
%box on

X = [-100;100;100;-100;-100];
Y = [-100;-100;100;100;-100];
Z = repmat(-100,5,1);
% draw a square in the xy plane with z = -100
plot3(X,Y,Z,'k','LineWidth',2);  
% draw a square in the xy plane with z = 100
plot3(X,Y,Z+200,'k','LineWidth',2); 
set(gca,'View',[-51.6 15.8]); % set the azimuth and elevation of the plot
for k=1:length(X)-1
   plot3([X(k);X(k)],[Y(k);Y(k)],[-100;100],'k','linewidth',2);
end

xlabel('x [$R_G$]','Interpreter','LaTex')
ylabel('y [$R_G$]','Interpreter','LaTex')
zlabel('z [$R_G$]','Interpreter','LaTex')
set(gca,'FontSize',18)
hold off