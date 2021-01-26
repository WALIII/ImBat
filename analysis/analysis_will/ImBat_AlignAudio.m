function audio_aligned = ImBat_AlignAudio(flightPaths,audio,cluster)


% 
dayIndex = find(contains(flightPaths.Dates,audio.session));


% Get all flights on this day and align the first cluster ( to do: all 3)

flightIndex2 = find(flightPaths.day == dayIndex & flightPaths.id == cluster );

for i = 1: length(flightIndex2);
flight_times(i) = flightPaths.AllFlightsTime(flightPaths.flight_starts_idx(flightIndex2(i)));
end
flight_times = sort(flight_times);
% now, find closest times:

disp('Finding closest times');
for i = 1:length(flight_times)
[minValueStart(i),closestIndexStart(i)] = min(abs(audio.time_vect -flight_times(i)));
end

% now, align audio:
clear RMS audio_aligned
counter = 1;
for ii = 1: length(closestIndexStart)-1;
 try
disp(['aligning audio to flight', num2str(ii), ' of ', num2str(length(closestIndexStart))]);
audio_aligned.audio_aligned(:,counter) = audio.data(closestIndexStart(ii)-(250000/audio.fs):closestIndexStart(ii)+1000000);
% RMS(:,ii) =zftftb_rms(audio_aligned(:,ii),audio.fs);
audio_aligned.RMS(:,counter) = downsample(tsmovavg(rms(abs(audio_aligned(:,ii)),2),'s',500,1),100);
counter = counter+1;
 catch
 disp('Flight is beyond the set pad');
end
end

audio_aligned.flight_times = flight_times;



