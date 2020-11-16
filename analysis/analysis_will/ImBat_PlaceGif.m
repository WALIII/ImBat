function ImBat_PlaceGif

axis tight
h = gcf;
filename = 'TestGif.gif';
grid on;
for i = 1:1:360
view(i,10)
axis tight
 % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif','DelayTime',0.05, 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append'); 
      end       


pause(0.1)
end