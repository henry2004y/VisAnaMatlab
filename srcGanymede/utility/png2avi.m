% script for creating avi
% Issue: this cannot deal with discontinous snapshots.
% Fix this later!
%
% Hongyang Zhou, hyzhou@umich.edu

writerObj = VideoWriter('../GroupMeeting/pressureRatio.avi');
writerObj.FrameRate=10;
open(writerObj);
for K = 1 : 50
  filename = sprintf('../GroupMeeting/PressureRatio/%4.4d.png', K);
  thisimage = imread(filename);
  writeVideo(writerObj, thisimage);
end
close(writerObj);