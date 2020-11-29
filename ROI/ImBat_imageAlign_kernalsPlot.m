figure('units','normalized','outerposition',[0 0 1 0.7]);
sgtitle('gaussfilt 70cutoff nofilt noexp b4 filt exp3 after');
filtKernal = [1 5 15 25 35 45 50];

IM1_aligned = imagesAligned.IM1_aligned;
IM2_aligned = imagesAligned.IM2_aligned;
IM3_aligned = imagesAligned.IM3_aligned;

for filt_i = 1:length(filtKernal)
hFilt = fspecial('disk',filtKernal(filt_i));
IM1_aligned_filt = imfilter(IM1_aligned,hFilt,'replicate');
IM2_aligned_filt = imfilter(IM2_aligned,hFilt,'replicate');
IM3_aligned_filt = imfilter(IM3_aligned,hFilt,'replicate');
IM1_aligned_filt = IM1_aligned_filt.^3;
IM2_aligned_filt = IM2_aligned_filt.^3;
IM3_aligned_filt = IM3_aligned_filt.^3;


%RGB the images
[aOverlap,bOverlap] = CaBMI_XMASS(IM1_aligned_filt,IM2_aligned_filt,IM3_aligned_filt,'hl',[.05 .3]);
%rgb the original unaligned for comparison
[aUnaligned,bUnaligned] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);

%quantify quartiles of overlap to see how good the days overlap
sizeIm = size(IM1_aligned_filt)/2;
quartRange(1,:) = [1,sizeIm(1),1,sizeIm(2)];
quartRange(2,:) = [1,sizeIm(1),sizeIm(2)+1,sizeIm(2)*2];
quartRange(3,:) = [sizeIm(1)+1,sizeIm(1)*2,1,sizeIm(2)];
quartRange(4,:) = [sizeIm(1)+1,sizeIm(1)*2,sizeIm(2)+1,sizeIm(2)*2];

for quart_i = 1:4
    %quartRange = [1+sizeIm(1)*(quart_i-1),quart_i*sizeIm(1),1+sizeIm(2)*(quart_i-1),quart_i*sizeIm(2)];
overlapQuart(quart_i,1) = sum(sum(abs(IM1_aligned_filt(quartRange(quart_i,1):quartRange(quart_i,2),quartRange(quart_i,3):quartRange(quart_i,4))-IM2_aligned_filt(quartRange(quart_i,1):quartRange(quart_i,2),quartRange(quart_i,3):quartRange(quart_i,4)))));
overlapQuart(quart_i,2) = sum(sum(abs(IM2_aligned_filt(quartRange(quart_i,1):quartRange(quart_i,2),quartRange(quart_i,3):quartRange(quart_i,4))-IM3_aligned_filt(quartRange(quart_i,1):quartRange(quart_i,2),quartRange(quart_i,3):quartRange(quart_i,4)))));
%overlapQuart(quart_i,3) = sum(sum(abs(IM1_aligned_filt(quartRange(quart_i,1):quartRange(quart_i,2),quartRange(quart_i,3):quartRange(quart_i,4))-IM3_aligned_filt(quartRange(quart_i,1):quartRange(quart_i,2),quartRange(quart_i,3):quartRange(quart_i,4)))));
medOverlapQuart(quart_i) = median(overlapQuart(quart_i,:));
end
medOverlapAll = median(overlapQuart,'all');

subplot(2,4,filt_i);
image(aOverlap);
set(gca,'xticklabel',[],'yticklabel',[]);
title(['rad ' num2str(filtKernal(filt_i)) ':' num2str(medOverlapAll)]);
end