function ImBat_PlaceCells(neuron, out);

% get general data on ROI traces

offset = 0.1; % account for slow calcium


% optional: add jitter to location:
for i = 1:3;
out.Location_jitter(:,i) = smooth(out.Location(:,i)+ randi(100,length(out.Location),1)/1,10);
out.Location(:,i) = smooth(out.Location(:,i),10);
end

% Smooth Neurons


% for i = 1:(size(out.Location,1)-1)
%     if out.Location(i,1) == out.Location(i+1,1) && out.Location(i,2) == out.Location(i+1,2);
%         for ii = 1:3
%         out.Location2(i,ii) = NaN;
%         end
%     else
%          for ii = 1:3
%         out.Location2(i,ii) = out.Location(i,ii);
%          end
%     end
% end
%         
    
% get 
figure(); 
hold on;
a1 = subplot(10,1,1); 
hold on;
% for i = 1:3
% plot(out.Location_time, out.Location(:,i)); 
% end


D = diff(nanmean(out.Location2(:,1:3),2));

plot(out.Location_time(1:end-1),(abs(D)))

a2 = subplot(10,1,2); 
hold on;
plot(out.video_times(1:end-1),zscore(smooth(mean((full(neuron.S(1:50,1:end-1))),1),100)));
ylim([0 5]);
a3 = subplot(10,1,3:10);
hold on;
for i = 1: 40; 
    hh = movmean(neuron.C_raw(i,:),4);
plot(out.video_times,zscore(hh)+i*4);
end
  linkaxes([a1, a2, a3], 'x');
 
  
  
%   figure(); % with heatmap
%   b1 = subplot(211); 
% hold on;
% for i = 1:3
% plot(out.Location_time, out.Location(:,i)); 
% end
% b2 = subplot(212); 
% hold on;
% h = imagesc(neuron.C);
% set(h, 'XData', out.video_times);
%   linkaxes([b1, b2], 'x');
  
  
  
  % plot an example trial
  
  % get min of location time
  [e e2] = min(out.Location_time);
  %
  
  
  
  % Plot 3d Scatter
  
  
  figure(); 

  for ii = 1:93; 
      hold on;
plot3(out.Location2(:,1),out.Location2(:,2),out.Location2(:,3),'k');
  
[~,xy] = find(neuron.S(ii,:)>0.1);  % get time
Spike_times = out.video_times(xy)-offset;

% find closeset location
% Spike_times

try
for i = 1:size(Spike_times,1)
try
    [minValue(:,i),closestIndex(:,i)] = min(abs(out.Location_time-Spike_times(i)));
LX(i) = out.Location2(closestIndex(:,i),1);
LY(i) = out.Location2(closestIndex(:,i),2);
LZ(i) = out.Location2(closestIndex(:,i),3);

catch
continue
end
end

disp([num2str(size(LX)),' Bursts in flight'])
scatter3(LX,LY,LZ,100,'or','filled');
title(['Cell no ',num2str(ii)]);
catch
    disp('cell not active');
    
    continue
end
   pause();
clf

clear LX LY LZ closestIndex Spike_times
  end
  
  
  
%   col = hsv(150);
%   
% figure(); 
% 
%   for ii = 1:40; 
%       hold on;
% plot3(out.Location(:,1),out.Location(:,2),out.Location(:,3),'k');
%   
% [~,xy] = find(neuron.S(ii,:)>0.1);  % get time
% Spike_times = out.video_times(xy);
% 
% % find closeset location
% % Spike_times
% 
% for i = 1:size(Spike_times,1)
% [minValue(:,i),closestIndex(:,i)] = min(abs(out.Location_time-Spike_times(i)));
% LX(i) = out.Location(closestIndex(:,i),1);
% LY(i) = out.Location(closestIndex(:,i),2);
% LZ(i) = out.Location(closestIndex(:,i),3);
% end
% 
% try
%     scatter3(LX,LY,LZ,20,col(ii,:));
% catch
%     disp('cell not active');
% end
% pause();
% clf
% 
% clear LX LY LZ closestIndex Spike_times
%   end
  
  