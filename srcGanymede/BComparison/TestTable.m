% Adding spacecraft location info underneath magnetic field plots.
% Written to be run after B_compare_steady.
%
% Hongyang Zhou, hyzhou@umich.edu 06/28/2018

%%
%format bank
% Create table for G8 B figure
data = {-1.55,-1.48,-1.37,-1.21,-1.05;-3.08,-1.14,0.83,2.77,4.70;...
   0.67,0.73,0.77,0.78,0.79};

for i = 1:numel(data)
   data{i} = sprintf('%3.2f', data{i});  
end


f = figure('Position', [100, 100, 1000, 700]);
h = subplot(511); LW1 = 2.5; LW2 = 1.5; FS=16; Height = 0.18;
plot(time, Bobs(:,1),'k',timesim,Bsim(:,1),'r','LineWidth',LW1); %Bx
h.XTickLabel = [];
h.Position(4) = Height;
ylabel('Bx [nT]');
legend(h,{'Obs','Sim'})
xlim([min(time) max(time)]);
h.XLim(1) = datetime(2000,5,20,9,45,0);
ylim([min(Bobs(:,1))-50 max(Bobs(:,1))+50 ]);
%title({'(a)',strcat('Galileo G',int2str(flyby),' Flyby Magnetic field')})
title(strcat('Galileo G',int2str(flyby),' Flyby Magnetic field'))
set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
h = subplot(512);
plot(time, Bobs(:,2),'k',timesim,Bsim(:,2),'r','LineWidth',LW1); %By
h.XTickLabel = [];
h.Position(4) = Height;
ylabel('By [nT]');
xlim([min(time) max(time)]);
h.XLim(1) = datetime(2000,5,20,9,45,0);
ylim([min(Bobs(:,2))-50 max(Bobs(:,2))+50 ]);
set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
h = subplot(513);
plot(time, Bobs(:,3),'k',timesim, Bsim(:,3),'r','LineWidth',LW1); %Bz
h.XTickLabel = [];
h.Position(4) = Height;
ylabel('Bz [nT]');
xlim([min(time) max(time)]);
h.XLim(1) = datetime(2000,5,20,9,45,0);
ylim([min(Bobs(:,3))-50 max(Bobs(:,3))+50 ]);
set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);
h = subplot(514);
plot(time,BobsStrength,'k',timesim,BsimStrength,'r','LineWidth',LW1); %|B|
h.Position(4) = Height;
ylabel('B [nT]');
xlim([min(time) max(time)]);
h.XLim(1) = datetime(2000,5,20,9,45,0);
ylim([min(BobsStrength)-50 max(BobsStrength)+50 ]);
set(gca,'FontSize',FS,'yMinorTick','on','LineWidth',LW2);


%f = figure;
t = uitable(f);
t.Data = data;
t.Position = [100,45,800,55];
t.ColumnName = '';
t.ColumnWidth = {150 150 150 150 150};
t.RowName = {'X','Y','Z'};


% Create textbox for G8
annotation(f,'textbox',[0.11 0.0915367020170514 0.812 0.07],...
   'String',{'X                -1.55                                                     -1.48                                                      -1.37                                                      -1.21                                                     -1.05','Y                -3.08                                                     -1.14                                                       0.83                                                        2.77                                                      4.70','Z                 0.67                                                       0.73                                                       0.77                                                        0.78                                                      0.79'},...
   'FitBoxToText','off',...
   'EdgeColor',[1 1 1]);


% Create textbox for G28
annotation(f,'textbox',[0.11 0.0915367020170514 0.812 0.07],...
   'String',{'X                       -0.71                                           -0.98                                         -1.23                                          -1.41                                         -1.57                                         -1.73','Y                       -5.21                                           -2.69                                         -0.14                                           2.44                                           4.94                                          7.46','Z                       -0.43                                           -0.43                                         -0.43                                          -0.40                                          -0.37                                        -0.33'},...
   'FitBoxToText','off',...
   'EdgeColor',[1 1 1]);