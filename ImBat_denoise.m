function FiltA =  ImBat_denoise(video);% % Video recon
% ImBat_denoise.m

% Remove wireless artifacts from CaIm video data

% WAL3
% d05/10/2019

display_mov = 0;
% Additional artifact rejection:
FiltA = double(video);

% remove initial baseline trend:
baseline = (smooth(squeeze(mean(mean(FiltA(1:20,:,:),1),2)),50));
baseline2 = baseline-mean(baseline);

for ii = 1: size(FiltA,3);
    FiltA(:,:,ii) =  FiltA(:,:,ii)-baseline2(ii);
end

figure();
for i = 1:10; % 100 itterations...
thresh = 3;
% remove large offsets

siga = squeeze(mean(mean(FiltA(:,1:20,:),1),2));
sigb = squeeze(mean(mean(FiltA(1:20,:,:),1),2));
sig = (siga+sigb)/2;
% sig_base = smooth(sig,50);
% sig_base = sig_base-mean(sig_base);
sig = zscore(sig);
sig = detrend(sig);
[k1 k2] = find(sig>thresh);
hold on;
plot(sig);
plot(k1,k2,'*')
disp(['Cleaning up ',num2str(size(k1,1)),' frames...']);
if size(k1,1) ==0;
    break
end

for ii = 1:size(k1,1)
    if k1(ii)<5

      FiltA(:,:,k1(ii)) = (median(FiltA(:,:,k1(ii)+2:k1(ii)+5),3));

    elseif k1(ii)>size(k1,1)-5
           FiltA(:,:,k1(ii)) = (median(FiltA(:,:,k1(ii)-4:k1(ii)-2),3));

    else

    FiltA(:,:,k1(ii)) = (median(FiltA(:,:,k1(ii)-4:k1(ii)-2),3)+median(FiltA(:,:,k1(ii)+2:k1(ii)+5),3))/2;

    end
end
end

% detrend data:

baseline = (smooth(squeeze(mean(mean(FiltA(1:20,:,:),1),2)),50));
baseline2 = baseline-mean(baseline);

for ii = 1: size(FiltA,3);
    FiltA(:,:,ii) =  FiltA(:,:,ii)-baseline2(ii);
end
figure(); plot(smooth(squeeze(mean(mean(FiltA(1:20,:,:),1),2)),50));

% % smooth video:
disp('smoothing video');
b = 3; %(smoothing )
% FiltB = mat2gray(FiltA-FBM);
FiltA = (convn(FiltA, single(reshape([1 1 1] / b, 1, 1, [])), 'same'));
mG = mean(video,3);
mF = mean(FiltA,3);


if display_mov ==1;
figure();
for i = 1:3:size(FiltA,3)
    colormap(gray)
 subplot(1,2,1)
imagesc(double(video(:,:,i))-mG,[-20 60]);

    subplot(1,2,2)
imagesc(FiltA(:,:,i)-mF,[-20 60]);
title('raw video');

%     subplot(1,3,3)
% imagesc(FiltB(:,:,i),[0.2 0.7]);
% title('raw video');

pause(0.01);
end
end
