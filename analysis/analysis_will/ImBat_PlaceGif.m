function ImBat_PlaceGif

% axis tight
axis off
h = gcf;
h.Color = 'w';
filename = 'TestGif_new.gif';
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
axis off

% go from view (1,1:10) first:
for i = flip(10:2:90)
view(1,i)
xlim([-2900 2900])
ylim([-2900 2900]);
% axis tight
 % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 90
          imwrite(imind,cm,filename,'gif','DelayTime',0.05, 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append'); 
      end       

pause(0.1)
end

for i = 1:1:360
view(i,10)
xlim([-2900 2900])
ylim([-2900 2900]);
% axis tight
 % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
%       if i == 1 
%           imwrite(imind,cm,filename,'gif','DelayTime',0.05, 'Loopcount',inf); 
%       else 
          imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append'); 
%       end       


pause(0.1)
end


% go from view (1,1:10) first:
for i = (10:2:90)
view(1,i)
xlim([-2900 2900])
ylim([-2900 2900]);
% axis tight
 % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
%       if i == 1 
%           imwrite(imind,cm,filename,'gif','DelayTime',0.05, 'Loopcount',inf); 
%       else 
          imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append'); 
%       end       

pause(0.1)
end