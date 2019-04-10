function timeline_method2(list)
  %make sure the list is ordered in reverse chronological order (to help identify the last label of a row)
  [~, order] = sortrows(vertcat(list{:, 2}), [-1 2]);
  list = list(order, :);
  %identify unique label and generate vertical positioning
  [labels, idx, ylabels] = unique(list(:, 1), 'stable');
  ypatches = max(ylabels) + 1 - [ylabels, ylabels, ylabels-1, ylabels-1]'; 
  ylabels = max(ylabels) + 1 - ylabels(idx);
  %generate horizonal positioning
  xpatches = [vertcat(list{:, 2}), fliplr(vertcat(list{:, 2}))]';
  xlabels = xpatches(2, idx);
  %plot
  figure;
  colour = parula(size(list, 1));
  %colour = parula(1);
  patch(xpatches, ypatches, reshape(colour, 1, [], 3));
  text(xlabels+5, ylabels+0.5, labels, 'fontsize', 10);
  xlabel('Time [s]');
  grid on
end

