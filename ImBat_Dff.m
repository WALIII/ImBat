function [Ymax, Y, maxFig] = ImBat_Dff(Y);
% ImBat_Dff

% Make df/f image

% Filter movie

Y = (convn(Y, single(reshape([1 1 1] /10, 1, 1, [])), 'same'));

% Take median of movie
Y_med = median(Y,3);

% Subtract median
Y = Y-Y_med;

% take max
Ymax = max(Y,[],3);
Ymax = imresize(Ymax,8);
maxFig = figure();
colormap(gray);
imagesc(Ymax);
hold on;







