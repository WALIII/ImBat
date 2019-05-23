function ImBat_PlaceCells(neuron, out);
% ImBat_PlaceCells

% Basic plotting of the location that certian cells are active...





% get general data on ROI traces

offset = 0.1; % account for slow calcium estimation ~move locations back 100ms in time... This is the knob to turn for 'prospective' coding...




% optional: add jitter to location:
for i = 1:3;
out.Location_jitter(:,i) = smooth(out.Location(:,i)+ randi(100,length(out.Location),1)/1,10);
out.Location(:,i) = smooth(out.Location(:,i),10);
end



%Using subplots, show:
% 1: the velocity/accel of the bat
% 2: average population activity ( smoothed, average of the
    % z-score-normalized data- will be usefull for finding 'replay' events
    % later on...
% 3: plotting individual neuron timeeries (z-scored) ( top 40)

figure(); 
hold on;

a1 = subplot(10,1,1); 
hold on;
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
 
  
  

  
 
% Plot the location in space that each cell is active

  mkdir('PlaceCells/fig'); % save as a figure file in local dir
   mkdir('PlaceCells/jpg'); % Save as .jpg in local dir
  figure(); 

  for ii = 1:93; % for each cell
      hold on;
plot3(out.Location2(:,1),out.Location2(:,2),out.Location2(:,3),'k');% plot the flight trajectory in space
  
[~,xy] = find(neuron.S(ii,:)>0.1);  % get time neuron is active
Spike_times = out.video_times(xy)-offset; % convert this to 'spike time'


try % this 'try/catch' is to avoid crashing if cells are not active in plotting window...
for i = 1:size(Spike_times,1)
try
    % Find the closest 'Location time' to the 'Spike time'
[minValue(:,i),closestIndex(:,i)] = min(abs(out.Location_time-Spike_times(i)));
LX(i) = out.Location2(closestIndex(:,i),1);
LY(i) = out.Location2(closestIndex(:,i),2);
LZ(i) = out.Location2(closestIndex(:,i),3);
catch % we need this if Spiketime occurs before/after the location tracking was on..
continue
end
end

% display, in the title, how many bursts there were:
disp([num2str(size(LX)),' Bursts in flight'])
scatter3(LX,LY,LZ,100,'or','filled');
title(['Cell no ',num2str(ii),'  ',num2str(size(LX)),' Bursts in flight']);
catch % if cell was not active...
    disp('cell not active');
    continue
end

% Save 'place cells' as jpg and fig files..
saveas(gcf,['PlaceCells/fig/','Cell_',num2str(ii)]);
saveas(gcf,['PlaceCells/jpg/','Cell_',num2str(ii),'.jpg']);

clf

% Clear the buffer for the next cell:
clear LX LY LZ closestIndex Spike_times
  end
  
  