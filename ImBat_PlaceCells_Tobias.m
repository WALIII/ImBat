function [plotFiringTrajectory] = ImBat_PlaceCells_Tobias(flightPaths, cellData, out)

global batName dateSesh sessionType topROI

offset = 0.1; % account for slow calcium estimation ~move locations back 100ms in time... This is the knob to turn for 'prospective' coding...
topROI = 30;

% Plot the location in space that each cell is active
plotFiringTrajectory =  figure();

for ii = 1:topROI; % for each cell
    hold on;
    plot3(out.flights(:,1),out.flights(:,2),out.flights(:,3),'k','LineWidth',2);% plot the flight trajectory in space
    
    [~,xy] = find(cellData.S(ii,:)>8);  % get time neuron is active
    Spike_times = out.video_times(xy)-offset; % convert this to 'spike time'
    
     
    try % this 'try/catch' is to avoid crashing if cells are not active in plotting window...
        for i = 1:size(Spike_times,1)
            try
                % Find the closest 'Location time' to the 'Spike time'
                [minValue(:,i),closestIndex(:,i)] = min(abs(out.Location_time-Spike_times(i)));
                LX(i) = out.flights(closestIndex(:,i),1);   
                LY(i) = out.flights(closestIndex(:,i),2);
                LZ(i) = out.flights(closestIndex(:,i),3);
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

