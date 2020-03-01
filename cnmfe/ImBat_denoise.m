function video =  ImBat_denoise(video,varargin);% % Video recon
% ImBat_denoise.m

% Remove wireless artifacts from CaIm video data

% WAL3
% d05/10/2019

% default params
display_mov = 0;
maxItter = 10;
thresh = 3;
smth_frames = 3;
ds_temp = 1; % temporal downsample ( default = 1, i.e. no downsampling);


% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'thresh'
            thresh=varargin{i+1};
        case 'itterations'
            maxItter = varargin{i+1};
        case 'smth_frames' % temporal downsample
            smth_frames = varargin{i+1};
        case 'metadata' % temporal downsample
            metadata = varargin{i+1};
            smth_frames =  metadata.initial_median_filter_kernal;
            smth_frames2 =  metadata.median_filter_kernal;
    end
end




% Additional artifact rejection:
video = single(video);

% remove initial baseline trend:
baseline = (smooth(squeeze(mean(mean(video(1:20,:,:),1),2)),50));
baseline2 = baseline-mean(baseline);

for ii = 1: size(video,3);
    video(:,:,ii) =  video(:,:,ii)-baseline2(ii);
end


% Median filter:
% % smooth video:
disp('Initial median filtering of video');
tic
[video] = ImBat_Filter(single(video),smth_frames);
toc

if metadata.artifact_reject ==1;
    figure();
    for i = 1:maxItter; % 100 itterations...
        
        % remove large offsets
        siga = squeeze(mean(mean(video(:,[1:20 (end-20):end],:),1),2));
        sigb = squeeze(mean(mean(video([1:20 (end-20):end],:,:),1),2));
        sig = (siga+sigb)/2;
        % sig_base = smooth(sig,50);
        % sig_base = sig_base-mean(sig_base);
        sig = zscore(sig);
        sig = detrend(sig);
        sig = abs(sig);
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
                
                video(:,:,k1(ii)) = (median(video(:,:,k1(ii)+2:k1(ii)+5),3));
                
            elseif k1(ii)>size(k1,1)-5
                video(:,:,k1(ii)) = (median(video(:,:,k1(ii)-4:k1(ii)-2),3));
                
            else
                
                video(:,:,k1(ii)) = (median(video(:,:,k1(ii)-4:k1(ii)-2),3)+median(video(:,:,k1(ii)+2:k1(ii)+5),3))/2;
                
            end
        end
    end
else
%
    
end

% detrend data:

baseline = (smooth(squeeze(mean(mean(video(1:20,:,:),1),2)),50));
baseline2 = baseline-mean(baseline);

for ii = 1: size(video,3);
    video(:,:,ii) =  video(:,:,ii)-baseline2(ii);
end

% % smooth video:
disp('Median filtering video');
[video] = ImBat_Filter(single(video),smth_frames2);

% Plot mean pixel values...
figure(); plot(smooth(squeeze(mean(mean(video(1:20,:,:),1),2)),50));
mG = mean(video,3);
mF = mean(video,3);





if display_mov ==1;
    figure();
    for i = 1:3:size(video,3)
        colormap(gray)
        subplot(1,2,1)
        imagesc(double(video(:,:,i))-mG,[-20 60]);
        
        subplot(1,2,2)
        imagesc(video(:,:,i)-mF,[-20 60]);
        title('raw video');
        
        
        pause(0.01);
    end
end
