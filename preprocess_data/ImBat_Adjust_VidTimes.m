function ROI_Data = ImBat_Adjust_VidTimes(ROI_Data);



for i = 1:size(ROI_Data,2);
    fs = ROI_Data{i}.ROIs.results.Fs;

Video_times = ROI_Data{i}.Alignment.out.video_times;
Video_times2 = Video_times(1:(30/fs):end);
% check ROI size:
lastFrame = size(ROI_Data{i}.ROIs.results.S,2);
ROI_Data{i}.Alignment.out.video_times2 = Video_times2(1:lastFrame);
end