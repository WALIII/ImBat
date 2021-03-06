function [Y3] = ImBat_FilterForMovCorr(Y)
% to remove 'spots' of dust on the CMOS


rad1 = 10;
thresh = -0.3;

% Take mean of image:
G = mean(Y,3);
% find bad spots
BW = G<thresh;
BW = single(BW);
h = fspecial('gaussian',5, 5);
BW2 = imfilter(BW,h);
BW2 = BW2>0.01;
BW = BW2;



% X = mask2poly(BW);
% X = boundary(X);
% [row,col] = size(G);
% 
%   scalF = 2.5;
%   scalF2 = scalF-1;
%   
% BW2 = poly2mask(X(:,1)*(scalF-mean(X(:,1)))*scalF2 ,X(:,2)*scalF -(mean(X(:,2))*scalF2),row,col);


% % expand beyond
% BW = imresize(BW,[e1]);
% BW = BW(4:end-3,4:end-3);


% Replicate video matrix:

Y2 = Y;
%Y3 = Y; %(backup)
% median filter the shit outta it
Y2 =  medfilt3(Y2 ,[9 9 1]);
disp('Done median filtering data');

% replace only these pixels with the smoothed image
for i = 1: size(Y2,3)
    temp = Y(:,:,i);
    temp2 = Y2(:,:,i);
    temp(BW) = temp2(BW);
    Y3(:,:,i) = temp;
end

figure(); 
subplot(1,2,1);
imagesc(mean(Y,3),[-.5 .5]);
title('before filtering');
subplot(1,2,2)
imagesc(mean(Y3,3),[-.5 .5]);
colorbar();
title('after filtering');
