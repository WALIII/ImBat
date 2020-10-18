function [out_video, lost_frames] = ImBat_DroppedFrames(video)

% detect dropped frames as a dramatic decrease in pixel variance in a frame;




sig = squeeze(mean(var(video([1:20 (end-20):end],[1:20 (end-20):end],:),1),2));
figure(); plot(sig);
title(' ID lost frames');
lost_frames = find(sig<median(sig)/4);

% loop through and remove low values
for i = 1:size(lost_frames)
    video(:,:,lost_frames(i)) =  median(video(:,:,lost_frames(i)-6:lost_frames(i)-1),3);
    disp(['repairing frame ', num2str(lost_frames(i))]);
end

out_video = video;
sig = squeeze(mean(var(out_video([1:20 (end-20):end],[1:20 (end-20):end],:),1),2));
figure(); plot(sig);
title(' repaired frames');