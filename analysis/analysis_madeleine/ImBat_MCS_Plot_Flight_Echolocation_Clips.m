% Sanity check to see what flights sound like 

% Load the audio you want to listen to
load('audioConCat_4.mat');
% Pad the audio and add the segment
padd = NaN((size(flightPaths34.trajectoriesSpline,2) - size(flightPaths34.MicrophoneVect,1))*1600,1);
audioConCat_padded = [padd;audioConCat];
sound(audioConCat_padded(audio_segment)*3,192000);

flight_number = [120,121,122];%[5,6,7,8];
pp = 500;
if length(flight_number) == 1
    FP_start = flightPaths34.flight_starts_idx_sorted(flight_number);
    FP_end = flightPaths34.flight_ends_idx_sorted(flight_number);
else
    FP_start = flightPaths34.flight_starts_idx_sorted(flight_number(1));
    FP_end = flightPaths34.flight_ends_idx_sorted(flight_number(end));
end
figure(); subplot(2,1,1); hold on;
plot3(flightPaths34.trajectoriesSpline(1,(FP_start-pp:FP_end+pp)),flightPaths34.trajectoriesSpline(2,(FP_start-pp:FP_end+pp)),flightPaths34.trajectoriesSpline(3,(FP_start-pp:FP_end+pp)));
for i=FP_start-pp:FP_end+pp
    if flightPaths34.EcholocationVect_padded(i) == 0.5
        plot3(flightPaths34.trajectoriesSpline(1,i),flightPaths34.trajectoriesSpline(2,i),flightPaths34.trajectoriesSpline(3,i),'*r');
    end
end
audio_segment = (FP_start-pp)*1600:(FP_end+pp)*1600;
audio_segment_ds = (FP_start-pp):(FP_end+pp);

subplot(2,1,2);  hold on;
plot(downsample(audioConCat_padded(audio_segment),1600));
plot(flightPaths34.EcholocationVect_padded(audio_segment_ds),'*r');
for i=1:length(flight_number)
    xline(flightPaths34.flight_starts_idx_sorted(flight_number(i))-FP_start+pp,'g','LineWidth',2); 
    xline(flightPaths34.flight_ends_idx_sorted(flight_number(i))-FP_start+pp,'r','LineWidth',2); 
end
sound(audioConCat_padded(audio_segment)*100,192000);