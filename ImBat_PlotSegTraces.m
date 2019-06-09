function ImBat_PlotSegTraces(neuron, out, out_traj)


% Run after:
% [out_traj] =  ImBat_SegTrajectories(out.Location,out.Location_time)



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
plot(out.Location_time(1:end-1),(abs(D)),'g')
plot(out.video_times(1:end-1),zscore(smooth(mean((full(neuron.S(1:50,1:end-1))),1),100)),'r');

 for i = 1: size(out_traj.ClusterIndex{1},2)
     to_plot = out_traj.flight_starts_times(out_traj.ClusterIndex{1}(i)):out_traj.flight_ends_times(out_traj.ClusterIndex{1}(i));
     plot(to_plot, ones(1,length(to_plot)),'LineWidth',5);
     plot(out_traj.flight_starts_times(out_traj.ClusterIndex{1}(i)), ones(1,length(out_traj.flight_starts_times(out_traj.ClusterIndex{1}(i)))),'r*')
 end
 

ylim([0 5]);

a3 = subplot(10,1,2:10);
hold on;
for i = 1: 40; 
    hh = movmean(neuron.C_raw(i,:),4);
plot(out.video_times,zscore(hh)+i*4);
end
  linkaxes([a1, a3], 'x');
  
  
  

 figure();
 hold on;
 for i = 1: size(out_traj.ClusterIndex{1},2)
     to_plot = out_traj.flight_starts_times(out_traj.ClusterIndex{1}(i)):out_traj.flight_ends_times(out_traj.ClusterIndex{1}(i));
     plot(to_plot, ones(1,length(to_plot)),'LineWidth',5);
     plot(out_traj.flight_starts_times(out_traj.ClusterIndex{1}(i)), ones(1,length(out_traj.flight_starts_times(out_traj.ClusterIndex{1}(i)))),'r*')
 end
  