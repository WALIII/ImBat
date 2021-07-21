function [out, metrics] = ImBat_alignTimeStamps(audio,video,TS,Markers);
% ImBat_alignTimeStamps

% Align the Freedomscope analog input with 'Cortex' software timestamps,
% outputs flight data [x,y,t] where 't' is timewarped to fit the FreedomScope data

% WAL3
% d05/10/2019



% core fucntions:
%            ImBat_formatTracking.m

%paramss
fs = audio.rate;% audio framerate
Tsfs = 120; % frame rate of the tracking/syncing software

DS_factor = round(fs/Tsfs); % should be '400'


% get tracking data
disp('extracting tracking data');
[Location, Location2] = ImBat_formatTracking(Markers); % formated location data
Rewards = TS(:,1);
TS = TS(:,2);

% TO DO: check which channel the audio is. Currenlty hard coded at ch2:
% Downample the audio data (which is the TTL input to the framegrabber in
% the imaging setup) to the rate of Cortex, 120Hz. The sample rate of the
% audio data is 48000Hz.
if size(audio.data,2)>1
audio_data = audio.data(:,2);
else
audio_data =audio.data;
end
audio_data = downsample(audio_data,DS_factor);
audio_tv = downsample(audio.times,DS_factor);

% Make time vector for TS
a = 1:length(TS);
a = a/Tsfs;
max(a);


TS_tv = a;
TS2 = TS;



audio_data(audio_data<0) = 0;
audio_z = zscore(audio_data);
TS_z = zscore(smooth(TS2,10));
try
[Apks,Alocs] = findpeaks(audio_z,'MinPeakProminence',4,'MinPeakDistance',60);
[Bpks,Blocs] = findpeaks(TS_z,'MinPeakProminence',1,'MinPeakDistance',6);
catch
Apks = 1;
Alocs = 1;
Bpks = 1;
Blocs = 1;
end


% make better sigfor vizualizations...
audio_infer = zeros(1,length(audio_z))';
audio_infer(Alocs) = 1;

TS_infer = zeros(1,length(TS_z))';
TS_infer(Blocs) = 1;

% offset:
offsetA = Alocs(1); % audio offset
offsetB = Blocs(1); % TS offset

% offset times:
audio_tv_offset = audio_tv - audio_tv(offsetA); % subtract this value
TS_tv_offset = TS_tv-TS_tv(offsetB);


% Get discrete reward signals
[~,rewardLocs] = findpeaks(Rewards,'MinPeakProminence',2);

% Align to first peak
figure();
hold on;
plot(audio_tv_offset,audio_z,'r');
plot(TS_tv_offset,TS_z,'b');


% % Make the first video timestamp "time zero"
TS_tv_offset = TS_tv-audio_tv_offset(1)-TS_tv(offsetB);
audio_tv_offset = audio_tv_offset-audio_tv_offset(1);

figure();
hold on;
plot(audio_tv_offset,audio_infer+.01,'r');
plot(TS_tv_offset,TS_infer,'b');
title('inferred offsets');



% explort the markers:
out.Location = Location;
out.Location_time = TS_tv_offset';
% out.markers_resampled % to the frame rate of the video
% out.markers_resampled_time
out.video_times = video.times;
out.Location2 = Location2;
out.RewardVector = Rewards;
out.RewardTime = rewardLocs;

% metrics for plotting later
metrics.audio_tv_offset = audio_tv_offset;
metrics.audio_infer = audio_infer;
metrics.TS_tv_offset = TS_tv_offset;
metrics.TS_infer = TS_infer;


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
