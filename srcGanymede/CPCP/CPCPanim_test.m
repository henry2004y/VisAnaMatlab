% CPCP scan movie plot.
%
% Require preprocessed CPCP data.

v = VideoWriter('CPCP vs time.avi');
v.FrameRate = 10;
v.open

% create new figure with specified size
hfig = figure(1);
set(hfig,'position', [10, 10, 800, 200]) 
colormap(jet);

npict = 600;

for ipict=1:npict
   
   plot(CPCPt,'-*');
   xline(ipict,'k','Current time','LabelOrientation','horizontal');
   
   frame = getframe(gcf);
   
   set(gca,'nextplot','replacechildren');
   
   clf;
   writeVideo(v,frame);
end

v.close
close(1)