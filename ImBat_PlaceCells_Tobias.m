function [plotFiringTrajectory] = ImBat_PlaceCells_Tobias(flightPaths, cellData, alignment)

global batName dateSesh sessionType topROI

offset = 0.1; % account for slow calcium estimation ~move locations back 100ms in time... This is the knob to turn for 'prospective' coding...
spikeThresh = 8; %threshold for eliminating noise from the S matrix
topROI = 30;

% Plot the location in space that each cell is active
plotFiringTrajectory =  figure();

for ii = 1:length(cellData.results.S(:,1)); % for each cell
    hold on;
    plot3(alignment.out.flights(:,1),alignment.out.flights(:,2),alignment.out.flights(:,3),'k','LineWidth',2);% plot the flight trajectory in space
    
    [~,xy] = find(cellData.results.S(ii,:)>spikeThresh);  % get time neuron is active
    xy = xy(xy<length(alignment.out.video_times));
    Spike_times = alignment.out.video_times(xy)-offset; % convert this to 'spike time'
    
     
    try % this 'try/catch' is to avoid crashing if cells are not active in plotting window...
        for i = 1:size(Spike_times,1)
            try
                % Find the closest 'Location time' to the 'Spike time'
                [minValue(:,i),closestIndex(:,i)] = min(abs(alignment.out.Location_time-Spike_times(i)));
                LX(i) = alignment.out.flights(closestIndex(:,i),1);   
                LY(i) = alignment.out.flights(closestIndex(:,i),2);
                LZ(i) = alignment.out.flights(closestIndex(:,i),3);
            catch % we need this if Spiketime occurs before/after the location tracking was on..
                continue
            end
        end
        
        % display, in the title, how many bursts there were:
        LX_s = LX;
        LX_s(isnan(LX)) = [];
        disp([num2str(size(LX_s)),' Bursts in flight'])
        hold on
        scatter3(LX,LY,LZ,100,'or','filled');
        title(['Cell no ',num2str(ii),'  ',num2str(size(LX)),' Bursts in flight']);
    catch % if cell was not active...
        disp('cell not active');
        continue
    end
    
    % Save 'place cells' as jpg and fig files..
    %saveas(gcf,['PlaceCells/fig/','Cell_',num2str(ii)]);
    %saveas(gcf,['PlaceCells/jpg/','Cell_',num2str(ii),'.jpg']);
    

    
    % Clear the buffer for the next cell:
    clear LX LY LZ closestIndex Spike_times
    
    clf
end

