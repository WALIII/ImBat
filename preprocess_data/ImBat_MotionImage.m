function [Im] = ImBat_MotionImage(Y);
% Take the fist and last thousand or so framse, compare dff image for
% overlap.

% WAL3
% d03/12/2020

Last = size(Y,3)-3000;
numFrames = 2000; % frames from begining to end..
Ydff = Y-mean(Y,3);
BB = squeeze(mean(Ydff(:,:,Last-numFrames:Last-1),3));
AA = squeeze(mean(Ydff(:,:,1:numFrames),3));

Im(:,:,1) = mat2gray(AA(10:end-10,10:end-10)); Im(:,:,3) = mat2gray(AA(10:end-10,10:end-10)); Im(:,:,2) = mat2gray(BB(10:end-10,10:end-10));

figure(); 
subplot(131);
imagesc(Im(:,:,1),[0 1]);
grid on;
title('First 1000 frames');

subplot(132);
imagesc(Im(:,:,2),[0 1]);
grid on;
title('Last 1000 frames');

subplot(1,3,3);
imagesc(imfilter(Im*1.4,1));
title('overlay');

grid on;


% figure(); imagesc(mat2gray(AA)-mat2gray(BB),[-2 2]); colormap(fireice);