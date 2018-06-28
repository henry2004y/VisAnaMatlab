
%%
figure('position', [1, 1, 800, 300]) % create new figure with specified size 
subplot_tight(2,3,1);
contourf(yq,zq,uxv,'Linestyle','none'); colorbar; colormap jet;
axis tight equal
%caxis([-100 100]);
title('Ux [km/s]');

subplot_tight(2,3,2);
contourf(yq,zq,uyv,'Linestyle','none'); colorbar; colormap jet;
axis tight equal
%caxis([-60 60]);
title('Uy [km/s]');

subplot_tight(2,3,3);
contourf(yq,zq,uzv,'Linestyle','none'); colorbar; colormap jet;
axis tight equal
%caxis([-80 80]);
title('Uz [km/s]');

subplot_tight(2,3,4);
contourf(yq,zq,jv,'Linestyle','none'); colorbar; colormap jet;
axis tight equal
caxis([0.1 0.9])
title('j')

hold on;
[centroid,boundaries,FTEcountJ] = find_FTE(jv,'VarThreshold',threshold_j,...
   'AreaThreshold',40);


for k = 1 : FTEcountJ
   % switch to grid axis
   thisBoundary = (boundaries{k}-1)*dy + [ymin zmin];
   plot(thisBoundary(:,1), thisBoundary(:,2), 'k', 'LineWidth', 1);
   % Put the "blob number" labels on the "boundaries" grayscale image.
   thiscentroid = centroid(k,:)*dy + [zmin ymin];
   text(thiscentroid(2), thiscentroid(1), num2str(k), ...
      'FontSize', 10, 'FontWeight', 'Bold')
end
hold off;

subplot_tight(2,3,5);
contourf(yq,zq,uv,'Linestyle','none'); colorbar; colormap jet;
axis tight equal
%caxis([20 140]);
title('U [km/s]')

subplot_tight(2,3,6);
contourf(yq,zq,pv,'Linestyle','none'); colorbar; colormap jet;
axis tight equal
caxis([2 12]);
title('p [nPa]')

% hold on;
% [centroid,boundaries,FTEcountP] = find_FTE(pv,'VarThreshold',threshold_p,...
%    'AreaThreshold',40);

%pthres(ipict) = graythresh((pv-min(pv(:))) ./ (max(pv(:))-min(pv(:))));

% for k = 1 : numel(boundaries)
%    % switch to grid axis
%    thisBoundary = (boundaries{k}-1)*dy + [ymin zmin];
%    plot(thisBoundary(:,1), thisBoundary(:,2), 'k', 'LineWidth', 1);
%    % Put the "blob number" labels on the "boundaries" grayscale image.
%    thiscentroid = centroid(k,:)*dy + [zmin ymin];
%    text(thiscentroid(2), thiscentroid(1), num2str(k), ...
%       'FontSize', 10, 'FontWeight', 'Bold')
% end
% hold off;
