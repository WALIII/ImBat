function [out] = ImBat_alignTimeStamps(audio,video,TS,Markers);
% ImBat_alignTimeStamps

% Align the Freedomscope analog input with 'Cortex' software timestamps,
% outputs flight data [x,y,t] where 't' is timewarped to fit the FreedomScope data

% WAL3
% d05/10/2019



%paramss
fs = audio.rate;% audio framerate
Tsfs = 120; % frame rate of the tracking/syncing software

DS_factor = round(fs/Tsfs); % should be '400'


% get tracking data
disp('extracting tracking data');
[Location] = ImBat_formatTracking(Markers); % formated location data

TS = TS(:,2);
audio_data =audio.data;
audio_data = downsample(audio_data,DS_factor);
audio_tv = downsample(audio.times,DS_factor);

% Make time vector for TS
a = 1:length(TS);
a = a/Tsfs;
max(a);


TS_tv = a;
TS2 = TS;



audio_z = zscore(smooth(abs(audio_data)))-min(zscore(smooth(abs(audio_data))));
TS_z = zscore(smooth(TS2,10));

[Apks,Alocs] = findpeaks(audio_z,'MinPeakProminence',4,'MinPeakDistance',60);
[Bpks,Blocs] = findpeaks(TS_z,'MinPeakProminence',1,'MinPeakDistance',6);

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


%ur = resample(u,3,2);

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



% find the offset frame value, and subtact it

% Make offset relative to the audio ( audio starts at 0)



% explort the markers:

out.Location = Location;
out.Location_time = TS_tv_offset';
% out.markers_resampled % to the frame rate of the video
% out.markers_resampled_time
out.video_times = video.times;
save('Location_data','out','-v7.3');% play the video

% Smooth Location data
for i = 1:3
    out.Location2(:,i) = movmean(out.Location(:,i),100);
end

% getting location data:
disp('finalizing cursor')
for i = 1:(size(out.Location,1))-2
    if out.Location2(i,1) == out.Location2(i+2,1) && out.Location2(i,2) == out.Location2(i+2,2);
        for ii = 1:3
        out.Location2(i,ii) = NaN;
        end
    else
         for ii = 1:3
        out.Location2(i,ii) = out.Location2(i,ii);
         end
    end
end


save('Alignment.mat','out');