function [echolocation_idxs,echolocation_vector_DS] = ImBat_MCS_extract_echolocations(audioConCat,out)

% Break audio into 2minute segments. Analyze each segment. Poop out the
% vector of peaks. NON-Downsampled.

clear echolocation_idxs
fs = 192000; ds_factor=10;
seconds = size(audioConCat,1)/fs;
two_minutes = ceil(seconds/120);
echolocation_idxs = [];

for i=1:two_minutes-1
    clear pks locs;
    if i==two_minutes-1
        short_seg = audioConCat(120*(i-1)*fs+1:end);
    else
        short_seg = audioConCat(120*(i-1)*fs+1:(120*(i-1)+120)*fs);
    end
    %figure(); hold on; plot(short_seg); title("Manually selected audio subsegment");
    
    % Remove infinite values
    if num2str(max(short_seg)) == '1'
        infinite_values = find(short_seg == max(short_seg));
        ii = NaN(size(short_seg,1),1);
        ii(infinite_values) = 0.5;
        %figure(); hold on; plot(short_seg); plot(ii,'*r');
        short_seg(infinite_values) = 0;
    end
    
    % Remove artifacts based on too-high power in high frequencies?
    
    % Extract peaks from NON-downsampled audio using bandpass filter and RMS
    % threshold.
    
    % Oilbird filtering ideas: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5451837/
    
    
%     % TEST W SHORTER SEG
%     shorter_seg = short_seg(6000000:10000000);
%     
%     % Filter the raw signal without bandpass and see how the different filters look:
%     % Smooth with a movmedian filter
%     kk=20;
%     MOV_mic_data = movmedian(shorter_seg,kk);
%     %ABS_MOV_mic_data = downsample(tsmovavg(rms(abs(shorter_seg),2),'s',500,1),100);
%     RMS_mic_data = tsmovavg(rms(abs(shorter_seg),2),'s',500,1); 
%     figure('name',"Raw Data"); hold on; x1=subplot(4,1,1); hold on; title("Raw Data"); plot(shorter_seg); x2=subplot(4,1,2); hold on; title("MovMedian Filter"); plot(MOV_mic_data); x3=subplot(4,1,3); hold on; title("RMS"); plot(RMS_mic_data); linkaxes([x1 x2],'x');
%     subplot(4,1,4); hold on; title("Periodogram")
%     % Make periodograms
%     periodogram(shorter_seg); hold off;
%     
%     
%     % Filter the raw signal WITH bandpass and see how the different filters look:
%     clear b a pxx w periodogram
%     [b,a]=butter(8,2*[40e3 75e3]./fs,'bandpass');
%     conv_mic_data=filter(b,a,shorter_seg);
%     % Smooth the bandpassed data with a movmedian filter
%     kk=20;
%     MOV_mic_data = movmedian(conv_mic_data,kk);
%     %ABS_MOV_mic_data_bp = downsample(tsmovavg(rms(abs(conv_mic_data),2),'s',500,1),100);
%     RMS_mic_data_bp = tsmovavg(rms(abs(conv_mic_data),2),'s',500,1);
%     figure('name',"Bandpass Filtered Data"); hold on; x1=subplot(4,1,1); hold on; title("BP Filtered Data"); plot(conv_mic_data); x2=subplot(4,1,2); hold on; title("MovMedian Filter"); plot(MOV_mic_data); x3=subplot(4,1,3); hold on; title("RMS"); plot(RMS_mic_data_bp); linkaxes([x1 x2],'x');
%     subplot(4,1,4); hold on; title("Periodogram")
%     % Make periodogram:
%     periodogram(conv_mic_data); hold off;
%     
%     
%     % Plot spectrograms of raw shorter_seg data
%     figure('name', 'Raw Data'); 
%     [IMAGE,F,T] = fb_pretty_sonogram(shorter_seg./abs(max(shorter_seg)),fs,'low',2.9,'zeropad',0);
%     colormap(hot)
%     imagesc(T,F,log(abs(IMAGE)+1e+2));set(gca,'YDir','normal');
%     ylabel('kHz')
%     xlabel('time (s)');
%     title('Spectrogram of manually selected audio subsegment');
%     
%      % Plot spectrograms of bandpass filtered shorter_seg data
%     figure('name', 'BandPass Filtered'); 
%     [IMAGE,F,T] = fb_pretty_sonogram(conv_mic_data./abs(max(conv_mic_data)),fs,'low',2.9,'zeropad',0);
%     colormap(hot)
%     imagesc(T,F,log(abs(IMAGE)+1e+2));set(gca,'YDir','normal');
%     ylabel('kHz')
%     xlabel('time (s)');
%     title('Spectrogram of BandPass Filtered Data');
%     
%     % From Markowitz using Ofer's software.
%     % Extract peaks from m_AM
%     [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,gravity_center, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]= sap_features(shorter_seg,192000);
%     [pks,locs] = findpeaks(abs(m_AM),'MinPeakProminence',0.4);
%     locs = locs*44;
%     pkvect = NaN(1,length(shorter_seg));
%     pkvect(locs) = 0.01;
%     figure(); hold on; 
%     subplot(4,1,1); hold on; plot(conv_mic_data); plot(pkvect,'*r');
%     subplot(4,1,2); hold on; plot(RMS_mic_data); plot(pkvect,'*r');
%     subplot(4,1,4); hold on; 
%     [IMAGE,F,T] = fb_pretty_sonogram(shorter_seg./abs(max(shorter_seg)),fs,'low',2.9,'zeropad',0);
%     colormap(hot)
%     imagesc(T,F,log(abs(IMAGE)+1e+2));set(gca,'YDir','normal');
%     ylabel('kHz')
%     xlabel('time (s)');
%     title('Spectrogram of manually selected audio subsegment');
%     sound(shorter_seg*15,fs);

    %%#-part spectrogram plot
%     figure(); % make spectrogram..
%     ax1 = subplot(3,1,1)
% 
%     [b,a]=ellip(5,.2,80,[100]/(fs/2),'high'); 
%     %[b,a]=ellip(5,.2,80,[3515]/(fs/2),'high'); % filter above 3515
%     [IMAGE,F,T] = fb_pretty_sonogram(filtfilt(b,a,short_seg./abs(max(short_seg))),fs,'low',2.9,'zeropad',0);
%     colormap(hot)
%     imagesc(T,F,log(abs(IMAGE)+1e+2));set(gca,'YDir','normal');
%     ylabel('kHz')
%     xlabel('time (s)');
%     title('Spectrogram of manually selected audio subsegment');
% 
%     % now plot RMS below
%     [b,a]=butter(8,2*[40e3 80e3]./fs,'bandpass');
%     convy=filter(b,a,short_seg);
%     temp1 = zscore(zftftb_rms(convy',fs))*500;
%     temp1(temp1<10) = 0; % get rid of filter artifacts..
%     hold on;
% 
%     ax2 = subplot(3,1,2);
%     plot((30:size(temp1,2))/fs,temp1(30:end),'b','LineWidth',3)
%     %linkaxes([ax1,ax2],'x'); axis tight
% 
%     % plot the OG data
%     ax3 = subplot(3,1,3); hold on; 
%     plot((30:size(temp1,2))/fs,short_seg(30:end)); plot((30:size(temp1,2))/fs,pkvect(30:end),'*r');
%     %plot(pk_vect,'*r'); 
%     linkaxes([ax1,ax2,ax3],'x'); axis tight
% 
%     %Listen to segment:
%     sound(short_seg*3, fs);
    
    % Actual filter process:
    [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,gravity_center, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]= sap_features(short_seg,192000);
    [pks,locs] = findpeaks(abs(m_AM),'MinPeakProminence',0.4);
    locs = locs*44;
    pkvect = NaN(1,length(short_seg));
    pkvect(locs) = 0.01;

    % Re-align the locations to be for the full vector
    locs = locs + (i-1)*120*fs;
    echolocation_idxs = [echolocation_idxs,locs'];
end

iii = NaN(length(audioConCat),1);
for i=1:length(echolocation_idxs)
    iii(echolocation_idxs(i)-800:echolocation_idxs(i)+800) = 0.5;
end
echolocation_vector_DS = downsample(iii,1600);

end