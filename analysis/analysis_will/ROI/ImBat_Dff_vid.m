function [Ymax, df] = ImBat_Dff_vid(Y,fileName);
% ImBat_Dff
vidLimits = [55 75];
% Make df/f image

% Filter movie
Y = (convn(Y, single(reshape([1 1 1] /10, 1, 1, [])), 'same'));

% Take median or mean of movie
Y_med = median(Y,3);
Y_mean = mean(Y,3);

% Subtract median or mean for df
YmedSub = Y-Y_med;
df = Y - Y_mean; 

% take max
Ymax = max(df,[],3);
Ymax = imresize(Ymax,1);
Ym = Ymax +100;
fig = figure();
colormap(gray);
imagesc(Ym);
axis off;
saveas(fig,[fileName '.tif']);

df = mat2gray(df)*100;
mov = figure();
v = VideoWriter([fileName '_' num2str(vidLimits) '.avi']);
open(v);    


for i=1:4:(length(df)/4)
    
    imagesc(mean(df(:,:,i:i+3),3),vidLimits);
    axis off;
        colormap(gray);
        frame = getframe(mov);
   writeVideo(v,frame);
   %disp(i);
    pause(0.001);
 
end
close(v);






