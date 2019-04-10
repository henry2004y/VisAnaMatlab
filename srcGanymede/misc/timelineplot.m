% Timeline plot

%%
data_with_init_time = [ 
       1, 10, 5, 3 ;
       3, 10, 3, 9 ;
       7, 10, 4, 8 ;
       12,10, 2, 2 ];

h = barh(data_with_init_time, 'stack');
set(h(1), 'facecolor', 'none', 'EdgeColor', 'none'); % disable the color of the first column (init time)
set(gca, 'YTickLabel', {'proc 1', 'proc 2', 'proc 3', 'proc 4'} ); % change the y axis tick to your name of the process
axis ij; % Put the first row at top

%%
dnv = datenum([2016 01 00]) + cumsum(eomday(2016, 1:12));       % X-Axis Dates
loc1 = datenum([2016 01 01]) + randi(364, 1, 3);                % Create Data
loc2 = datenum([2016 01 01]) + randi(364, 1, 5);                % Create Data
loc3 = datenum([2016 01 01]) + randi(364, 1, 7);                % Create Data
figure(1)
plot(loc1, 3*ones(size(loc1)), 'x')
hold on
plot(loc2, 2*ones(size(loc2)), 'x')
plot(loc3, 1*ones(size(loc3)), 'x')
hold off
grid
set(gca, 'XTick', dnv)
datetick('x', 'mmm', 'keepticks')
axis([xlim    0  4])

%% use timeline_method2
%https://www.mathworks.com/matlabcentral/answers/264668-creation-of-timeline-with-multiple-occurrences-of-various-events
timeline2({'A', [-10 20]
           'B', [0 50]
           'C', [-100 100]
           'A', [90 100]
           'D', [30 50]
           'D', [-20 20]})
        
%%
timecount(timecount~=0) = 1;

iInterval = 1;
% J loop
for ipict = 2:length(timecount)-1
   if timecount(ipict,1) == 1
      if timecount(ipict-1,1) == 0
         first = ipict;
      elseif timecount(ipict+1,1) == 0
         last = ipict;
         J(iInterval,:) = [first last];
         iInterval = iInterval + 1;
      end
   else
      continue
   end  
end

iInterval = 1;
% Pe loop
for ipict = 2:length(timecount)-1
   if timecount(ipict,2) == 1
      if timecount(ipict-1,2) == 0
         first = ipict;
      elseif timecount(ipict+1,2) == 0
         last = ipict;
         Pe(iInterval,:) = [first last];
         iInterval = iInterval + 1;
      end
   else
      continue
   end
end

% timeline_method2({'J',[3  52],
%            'J',[68 72],
%            'J',[74 79],
%            'J',[91 94],
%            'J',[111 127],
%            'J',[239 243],
%            'J',[291 305],
%            'J',[539 543],
%            'J',[548 551],
%            'Pe',[11 52],
%            'Pe',[79 83],
%            'Pe',[95 100],
%            'Pe',[107 109],
%            'Pe',[116 125],
%            'Pe',[129 138],
%            'Pe',[152 154],
%            'Pe',[167 172],
%            'Pe',[204 206],
%            'Pe',[295 311],
%            'Pe',[316 326],
%            'Pe',[335 336],
%            'Pe',[338 339],
%            'Pe',[397 399],
%            'Pe',[402 405],
%            'Pe',[420 450],
%            'Pe',[452 479],
%            'Pe',[500 501],
%            'Pe',[544 553]})
              
%% Count by hand
timeline_method2({'G8',[0 3]
   'G8',[3 51]
   'G8',[62 65]
   'G8',[66 72]
   'G8',[73 83]
   'G8',[88 103]
   'G8',[106 140]
   'G8',[183 191]
   'G8',[204 209]
   'G8',[232 246]
   'G8',[273 277]
   'G8',[280 286]
   'G8',[287 314]
   'G8',[313 326]
   'G8',[352 356]
   'G8',[356 361]
   'G8',[376 384]
   'G8',[381 385]
   'G8',[409 414]
   'G8',[417 420]
   'G8',[497 502]
   'G8',[506 510]
   'G8',[512 521]
   'G8',[529 554]
   'G8',[560 566]
   'G8',[592 596]
   'G8',[599 600]
   'G28',[4 12]
   'G28',[16 34]
   'G28',[38 88]
   'G28',[93 113]
   'G28',[119 126]
   'G28',[128 132]
   'G28',[163 170]
   'G28',[198 220]
   'G28',[239 260]
   'G28',[293 304]
   'G28',[320 329]
   'G28',[334 366]
   'G28',[368 400]
   'G28',[426 430]
   'G28',[433 454]
   'G28',[460 479]
   'G28',[486 506]
   'G28',[530 538]
   'G28',[564 582]
   'G28',[587 597]})


%%
lineNames={'G8' 'G28'};
startTimes = {[0 3 62 66 73 88 106 183 204 232 273 280 287 313 352 356 376 381 409 417 436 464 486 497 506 512 529 560 592 599],...
   [4 16 37 93 119 128 163 198 239 293 320 334 426 433 460 486 530 564 587]};
endTimes = {[3 51 65 72 83 103 140 191 209 246 277 286 314 326 356 361 384 385 414 420 453 478 488 502 510 521 554 566 596 600],...
   [12 34 88 113 126 132 170 220 260 304 329 366 400 430 454 479 506 538 582 597]};
% Call timeline.m.
patchHndls = timeline(lineNames,startTimes,endTimes,'lineSpacing',.1,'facecolor','b');
%datetick('keeplimits'); 
title('FTE timeline chart');
xlabel('time [s]');
set(gcf,'position',[300 300 1000 200]);
set(gca,'FontSize',14,'LineWidth',1.2);
