function [out, metrics] = ImBat_MCS_alignTimeStamps(audio,video,TS,Markers,audioConCat,ttlConCat);
% ImBat_alignTimeStamps

% Align the Freedomscope analog input with 'Cortex' software timestamps,
% outputs flight data [x,y,t] where 't' is timewarped to fit the FreedomScope data

% WAL3
% d05/10/2019

% core fucntions:
%            ImBat_formatTracking.m

%paramss
thresh=0.0013;%voltage threshold for echolocation clicks for earthwork mics
thresh2=0.005;%voltage threshold for echolocation clicks for knowles mics
calldist=0.015;%minimal distance between two clicks (in seconds)
errorval=0.005;%added distance in seconds to allow for slight peak detection errors
fs_AV = audio.rate;% audio framerate
Tsfs = 120; % frame rate of the tracking/syncing software CORTEX

DS_factor_AV = round(fs_AV/Tsfs); % should be '400'

% get tracking data
disp('extracting tracking data');
[Location, Location2] = ImBat_formatTracking(Markers); % formated location data
Rewards = TS(:,1);
TS = TS(:,2);


% TO DO: check which channel the audio is. Currenlty hard coded at ch2:
if size(audio.data,2)>1
    audio_data = audio.data(:,2);
else
    audio_data = audio.data;
end
audio_data = downsample(audio_data,DS_factor_AV);
audio_tv = downsample(audio.times,DS_factor_AV);

fs_mic = 192000; % microphone frame rate
Tsfs = 120; % frame rate of the tracking/syncing software
DS_factor_mic = round(fs_mic/Tsfs); 
ttl_data_DS = downsample(ttlConCat,DS_factor_mic);
mic_data_DS = downsample(audioConCat,DS_factor_mic);

% LOOK IF lenght(TTL_DATA_DS) is different than audio_data by more than 12
% seconds)

seconds_of_mic_data = round(length(audioConCat)/192000);
seconds_of_cortex_data = round(size(Location,1)/120);

% Fill in the echolocation peak vector and then downsample to match others
% INSTEAD: find the time zero for the microphone timing vector THEN extract
% the echolocations.
% ii=NaN(size(audioConCat,1),1);
% for i=1:length(echolocation_peaks)
%     ii(echolocation_peaks(i)-800:echolocation_peaks(i)+800)=0.5;
% end
% echo_data_DS = downsample(ii,DS_factor_mic);

% % Plot to check
% st = 1000; ed = 2000;
% figure(); 
% subplot(2,1,1); hold on;
% plot(audioConCat(st*DS_factor_mic:ed*DS_factor_mic));
% iii=NaN(size(audioConCat,1),1);
% iii(echolocation_peaks) = 0.5;
% plot(iii(st*DS_factor_mic:ed*DS_factor_mic),'*r');
% subplot(2,1,2); hold on;
% plot(mic_data_DS(st:ed));
% plot(echo_data_DS(st:ed),'*r');

% Plot the beginning and end of the Mic TTL data to see if you can align to
% start or end.
figure(); 
subplot(2,1,1); hold on; 
title('Start of downsampled TTL');
plot(ttl_data_DS(1:1000))
subplot(2,1,2); hold on; 
title('End of downsampled TTL');
plot(ttl_data_DS(end-1500:end));

% Make time vector for TS (Cortex TTL)
a = 1:length(TS);
a = a/Tsfs;
max(a);

% Make time vector for ttl_data % Downsampled microphone TTL 
c = 1:length(ttl_data_DS);
c = c/Tsfs;
max(c);

TS_tv = a; 
ttl_data_tv = c; 
TS2 = TS;

audio_data(audio_data<0) = 0; ttl_data_DS(ttl_data_DS < 0) = 0; 
audio_z = zscore(audio_data); ttl_data_z = zscore(ttl_data_DS); 
TS_z = zscore(smooth(TS2,10)); 
try
[Apks,Alocs] = findpeaks(audio_z,'MinPeakProminence',4,'MinPeakDistance',60);
[Bpks,Blocs] = findpeaks(TS_z,'MinPeakProminence',1,'MinPeakDistance',6);
[Cpks,Clocs] = findpeaks(ttl_data_z,'MinPeakProminence',1,'MinPeakDistance',60);
catch
Apks = 1;
Alocs = 1;
Bpks = 1;
Blocs = 1;
Cpks = 1;
Clocs = 1;
end

%% START alignment case
% make better sigfor vizualizations... (analog signal has peaks at
% different hights, make it binary instead)
audio_infer = zeros(1,length(audio_z))';
audio_infer(Alocs) = 1;

TS_infer = zeros(1,length(TS_z))';
TS_infer(Blocs) = 1;

ttl_data_infer = zeros(1,length(ttl_data_z))';
ttl_data_infer(Clocs) = 1;

% offset:
offsetA = Alocs(1); % audio offset
offsetB = Blocs(1); % TS offset
offsetC = Clocs(1); % mic offset

% offset times
audio_tv_offset = audio_tv - audio_tv(offsetA); % subtract this value
TS_tv_offset = TS_tv-TS_tv(offsetB);
ttl_tv_offset = ttl_data_tv - ttl_data_tv(offsetC);

% % Align to first peak
% figure();
% hold on;
% plot(audio_tv_offset,audio_z,'r');
% plot(TS_tv_offset,TS_z,'b');
% plot(ttl_tv_offset,ttl_data_z,'g');

% offset times for aligning to START
audio_tv_offset = audio_tv - audio_tv(offsetA); % Originalf framegrabber TTL stream, subtract the index of the first peak.
TS_tv_offset = TS_tv - TS_tv(offsetB) - audio_tv_offset(1) - audio_tv(offsetA);% - audio_tv(offsetA); % Origianl Cortex TTL stream, subtract the index of the first peak, then the first index of the framegrabber TTL, then the first index of the framegrabber audio stream
ttl_tv_offset = ttl_data_tv - ttl_data_tv(offsetC) - audio_tv_offset(1) - audio_tv(offsetA);% - audio_tv(offsetA); %- audio_tv(offsetB);

% inferred offsets alining to start
figure();
hold on;
plot(audio_tv_offset,audio_infer+.01,'r');
plot(TS_tv_offset,TS_infer,'b');
plot(ttl_tv_offset,ttl_data_infer-.02,'g');
title('Inferred offsets aligned to START. Framegrabber red, Cortex blue, mic green.');

%% END Case
% offset for ALIGNING TO THE LAST TTL IN SPECIAL CASES
E_offsetA = Alocs(end); % audio offset
E_offsetB = Blocs(end); % TS offset
E_offsetC = Clocs(end); % mic offset

% offset times for aligning to END
E_audio_tv_offset = audio_tv - audio_tv(E_offsetA); % Originalf framegrabber TTL stream, subtract the index of the first peak.
E_TS_tv_offset = TS_tv - TS_tv(E_offsetB) - E_audio_tv_offset(1) - audio_tv(offsetA); % Origianl Cortex TTL stream, subtract the index of the first peak, then the first index of the framegrabber TTL, then the first index of the framegrabber audio stream
E_ttl_tv_offset = ttl_data_tv - ttl_data_tv(E_offsetC) - E_audio_tv_offset(1) - audio_tv(offsetA); %- audio_tv(offsetB);

E_audio_tv_offset = audio_tv - audio_tv(E_offsetA) - E_audio_tv_offset(1) - audio_tv(offsetA);

% aligning to END
figure();
hold on;
plot(E_audio_tv_offset,audio_infer+.01,'r');
plot(E_TS_tv_offset,TS_infer,'b');
plot(E_ttl_tv_offset,ttl_data_infer-.02,'g');
title('Inferred offsets aligned to END. Framegrabber red, Cortex blue, mic green.');

% Get discrete reward signals
[~,rewardLocs] = findpeaks(Rewards,'MinPeakProminence',2);

% Plot the distribution of times between TTLs to see if there are TTLs
% missing in the middle of the session or just at the end.
ttl_timing = diff(ttlConCat);
ttl_timing_infer = diff(ttl_data_infer);
ttl_timing_idx = find(ttl_timing_infer == 1);
ttl_timing_int = [ttl_timing_idx(1);diff(ttl_timing_idx)];
figure(); hold on; plot(ttl_timing_int,'o','MarkerFaceColor','r');
title("Time between TTLs. Range should be VERY small");
if range(ttl_timing_int) > 5
    high_o = find(ttl_timing_int > 365);
    low_o = find(ttl_timing_int < 355);
end
if ~isempty(high_o)
    if high_o == 1
        disp("The first TTL was after the start of the recording, good. Align to START");
        align_ = 'start';
    else
        disp("Extra long gap between TTLS: There is a TTL missing from in the middle");
        align_ = [];
    end
end
if ~isempty(low_o)
    if low_o==1
        disp("The first TTL was after the start of the recording, good. Align to END");
        align_ = 'end';
    else
        disp("Extra short gap between TTLS: There is a TTL missing from in the middle");
        align_ = [];
    end
end

if isempty(align_)
    disp("Look at TTLS! There is a mistake in the middle, perhaps catastophic data loss.");
elseif strcmp(align_,'start')
    disp("Aligning to START");
    out.Location_time = TS_tv_offset'; % WAS ORIGINALLY TS_tv_offset. MCS
    out.Microphone_Time = ttl_tv_offset'; % WAS ORIGINALLY ttl_tv_offset. MCS
elseif strcmp(align_,'end')
    disp("Aligning to END");
    out.Location_time = E_TS_tv_offset'; % WAS ORIGINALLY TS_tv_offset. MCS
    out.Microphone_Time = E_ttl_tv_offset'; % WAS ORIGINALLY ttl_tv_offset. MCS
end
%% Export the markers and everything in the out struct.
% Decide how to do the alignment
% export the markers:
out.Location = Location;
% out.markers_resampled % to the frame rate of the video
% out.markers_resampled_time
out.video_times = video.times;
out.Location2 = Location2;
out.RewardVector = Rewards;
out.RewardTime = rewardLocs;
out.MicrophoneVector = mic_data_DS;

% metrics for plotting later
metrics.audio_tv_offset = audio_tv_offset;
metrics.audio_infer = audio_infer;
metrics.TS_tv_offset = TS_tv_offset;
metrics.TS_infer = TS_infer;
metrics.ttl_infer = ttl_data_infer;
metrics.ttl_tv_offset = ttl_tv_offset';
metrics.E_audio_tv_offset = E_audio_tv_offset;
metrics.E_TS_tv_offset = E_TS_tv_offset';
metrics.E_ttl_tv_offset = E_ttl_tv_offset';

% Smooth Location data
 for i = 1:3
     out.Location3(:,i) = movmean(out.Location(:,i),30);
 end

out.flights = out.Location3; % take out the rest data...
% getting location data:
disp('finalizing cursor')
nancard = 0;
for i = 1:(size(out.Location,1))-20
    if out.Location3(i,1) == out.Location3(i+20,1) && out.Location3(i,2) == out.Location3(i+20,2);
        for ii = 1:3
            if nancard ==0;
                lastval = out.flights(i,ii);
                nancard =1;
            end
        out.flights(i,ii) = NaN;
        end
    else
         for ii = 1:3
             if nancard == 1;
                 out.flights(i,ii) = lastval;
                 nancard = 0;
             else
        out.flights(i,ii) = out.flights(i,ii);
             end
         end
    end
end


save('Alignment.mat','out','metrics','-v7.3');
