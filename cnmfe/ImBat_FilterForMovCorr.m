function [Y3] = ImBat_FilterForMovCorr(Y)
% to remove 'spots' of dust on the CMOS



thresh = -0.4;

% Take mean of image:
G = mean(Y,3);

% find bad spots
BW = G<thresh;

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
