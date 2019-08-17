function [Ymax, Y, maxFig] = ImBat_Dff(Y);

global batName dateSesh sessionType;

% ImBat_Dff
scaling = 8;
% Make df/f image

% Filter movie

Y = (convn(Y, single(reshape([1 1 1] /10, 1, 1, [])), 'same'));

% Take median of movie
Y_med = median(Y,3);

% Subtract median
Y = Y-Y_med;

% take max
Ymax = max(Y,[],3);
Ymax = imresize(Ymax,scaling);
maxFig = figure();
colormap(gray);
imagesc(Ymax);
hold on;
xticks = get(gca,'xtick');
yticks = get(gca,'ytick');
scaling  = 1.1; %1.1um per pixel
newlabelsX = arrayfun(@(ax) sprintf('%g', scaling * ax), xticks, 'un', 0);
newlabelsY = arrayfun(@(ay) sprintf('%g', scaling * ay), yticks, 'un', 0);
set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
title(['Max Projection dFF: ' batName ' ' dateSesh ' ' sessionType]);
xlabel('um'); ylabel('um');







