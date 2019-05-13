function ImBat_PlaceCells(neuron, out);

% get general data on ROI traces

figure(); 
hold on;
a1 = subplot(10,1,1:2); 
hold on;
for i = 1:3
plot(out.Location_time, out.Location(:,i)); 
end

a2 = subplot(10,1,3:4); 
hold on;
plot(out.video_times(1:end-1),smooth(mean((neuron.S(1:50,1:end-1)),1),100));

a3 = subplot(10,1,5:10);
hold on;
for i = 1: 40; 
plot(out.video_times,zscore(neuron.C_raw(i,:))+i*3);
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
  [e e2] = min(out.Location_time)
  %
  
  
  
  % Plot 3d Scatter
  
  col = hsv(10);
  
figure(); 
hold on;
plot3(out.Location(:,1),out.Location(:,2),out.Location(:,3),'k');
  
  for ii = 1:10 
[~,xy] = find(neuron.S(ii,:)>0.1);  % get time
Spike_times = out.video_times(xy);

% find closeset location
% Spike_times

for i = 1:size(Spike_times,1)
[minValue(:,i),closestIndex(:,i)] = min(abs(out.Location_time-Spike_times(i)));
LX(i) = out.Location(closestIndex(:,i),1);
LY(i) = out.Location(closestIndex(:,i),2);
LZ(i) = out.Location(closestIndex(:,i),3);
end

scatter3(LX,LY,LZ,20,col(ii,:));
clear LX LY LZ closestIndex Spike_times
  end
  
  