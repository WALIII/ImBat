function make_bat_task_video
%cd('/Volumes/tron/data/bat/newchill/190416/video')
d1 = dir('1_*');
d2 = dir('2_*');
for allv = 1 : size(d1,1)
   vid1 = VideoReader(d1(allv).name);
   vid2 = VideoReader(d2(allv).name);
   if allv == 1
       videoPlayer = vision.VideoPlayer;
       % new video
       outputVideo = VideoWriter('combined.avi');
       outputVideo.FrameRate = vid1.FrameRate;
       open(outputVideo);
   end
   while hasFrame(vid1) && hasFrame(vid2)
       img1 = readFrame(vid1);
       img2 = readFrame(vid2);
       imgt = horzcat(img1, img2);
       % play video
       step(videoPlayer, imgt);
       % record new video
       writeVideo(outputVideo, imgt);
   end
end
release(videoPlayer);
close(outputVideo);