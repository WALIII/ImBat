
function ImBat_PlotAllPix(FlightAlignedY);
% Plot all pix data for entire clustered series

filt2use = 10;
exp2use = 2;
upsamp_fact = 10;
numClust = size(FlightAlignedY,2);
h=fspecial('gaussian',10,10);
counter = 1;
for i = 1:numClust
    try
    YIM = FlightAlignedY{i};
    % filter and subtract background
    YIM1 = YIM;
    for ii = 1:size(YIM,4);
       %YIM1(:,:,:,ii) = medfilt3(squeeze(YIM(:,:,:,ii)),[1 1 1]);
       BG = min(squeeze(YIM1(:,:,:,ii)),[],3);
       BG=imfilter(BG,h,'circular');
       YIM1(:,:,:,ii) = squeeze(YIM1(:,:,:,ii))./(BG+.01);
       YIM1(:,:,:,ii) = medfilt3(squeeze(YIM(:,:,:,ii)),[3 3 3]);
       % Better version?
%         BG = min(squeeze(YIM1(:,:,:,ii)),[],3);
%        BG2 = max(squeeze(YIM1(:,:,1:125,ii)),[],3);
%        BG=imfilter(BG,h,'circular');
%        YIM1(:,:,:,ii) = squeeze(YIM1(:,:,:,ii))-BG2;
%        YIM1(:,:,:,ii) = (squeeze(YIM1(:,:,:,ii))./(BG+.01))*100;
%        YIM1(:,:,:,ii) = medfilt3(squeeze(YIM1(:,:,:,ii)),[3 3 3]);

    end
    
    % figure(); plot(mean(squeeze(mean(YIM_4PNR(1:40,1:40,:),1)))')
    YIM_4PNR = squeeze(median(YIM1(:,:,50:end-100),4));
    YIM1 = squeeze(median(YIM1,4));
    YIM2 = YIM1(:,:,100:450);
    
    [Cn, PNR, PNR_mov] = ImBat_correlation_image(imresize(YIM_4PNR,2));
    YIM_max = imresize(PNR,0.5);%std(YIM2,[],3);
    [im1_rgb norm_max_proj{counter},I{counter},idx_img{counter}] = CABMI_allpxs(imresize(YIM2,upsamp_fact),'filt_rad',filt2use,'exp',exp2use,'mask2use',imresize(double(YIM_max),upsamp_fact),'per',0);
    counter = counter+1;
    clear YIM YIM1
    catch
        disp('missing day, moving on');
    end
end

figure();
for i = 1:length(I);
    subplot(1,counter,i);
    imagesc(I{i});
end