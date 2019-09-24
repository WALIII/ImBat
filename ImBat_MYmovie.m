
function ImBat_MYmovie(out,Y)

%Y = imresize(Y,0.5)
% find the first timestamp of the location data

[a b] = min(find(out.Location(:,1)~=0))
t = out.Location_time(a);

% clean up movie:



[minValue,closestIndex] = min(abs(t-out.video_times));


fig = figure();
counter = 1;
mY = mean(Y(:,:,closestIndex:size(Y,3)),3);
Y = Y-mY;
% Y = mat2gray(Y);

v = VideoWriter('peaks6.avi');
v.FrameRate = 90;
open(v);

for i = closestIndex-300:closestIndex+4000;%size(Y,3);
    
    % step forward in video, find closest time in tracking
    [minValue,closestIndex2] = min(abs(out.video_times(i)-out.Location_time));
    
    fC(counter,:) = out.Location(closestIndex2,:);
    
    subplot(6,1,1:3)
    imagesc(imresize(Y(:,:,i),4),[-5 20]);
    colormap(gray);
    subplot(6,1,4:6);
    plot3(fC(:,1),fC(:,2),fC(:,3));
    xlim([-2600 2600]);
    ylim([-2600 2600]);
    zlim([-2600 2600]);
    grid on
    if counter ==1;
        disp('resize window!')
        pause(10);
    end
        counter = counter+1;

    frame = getframe(fig);
    
    writeVideo(v,frame);
    pause(0.001);
end

close(v);
