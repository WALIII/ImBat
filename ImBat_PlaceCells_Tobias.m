function [plotFiringTrajectory] = ImBat_PlaceCells_Tobias(flightPaths, cellData, alignment,varargin)

global topROI

offset = 0.1; % account for slow calcium estimation ~move locations back 100ms in time... This is the knob to turn for 'prospective' coding...
spikeThresh = 8; %threshold for eliminating noise from the S matrix
topROI = 30;

batName = [];
dateSesh = [];
sessionType = [];
saveFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'batname'
            batName=varargin{i+1};
        case 'datesesh'
            dateSesh = varargin{i+1};
        case 'sessiontype'
            sessionType = varargin{i+1};
        case 'saveflag'
            saveFlag = varargin{i+1};
    end
end

%labels for loading and saving data if running independent fromImBat_analyze
if saveFlag == 1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
    
    cellData = load([pwd '/processed/Motion_corrected_Data_DS_results.mat']);
    alignment = load([pwd '/processed/Alignment.mat']);
    load([pwd '/analysis/' label '_flightPaths.mat']);
end

% Plot the location in space that each cell is active
plotFiringTrajectory =  figure();

for ii = 1:length(cellData.results.S(:,1)); % for each cell
    hold on;
    plot3(alignment.out.flights(:,1),alignment.out.flights(:,2),alignment.out.flights(:,3),'k','LineWidth',2);% plot the flight trajectory in space
    
    
    [~,xy] = find(cellData.results.S(ii,:)>spikeThresh);  % get time neuron is active
    xy = xy(xy<length(alignment.out.video_times));
    Spike_times = alignment.out.video_times(xy)-offset; % convert this to 'spike time'
    peak_heights = cellData.results.S(ii,xy);
    
    
    try % this 'try/catch' is to avoid crashing if cells are not active in plotting window...
        for i = 1:size(Spike_times,1)
            try
                % Find the closest 'Location time' to the 'Spike time'
                [minValue(:,i),closestIndex(:,i)] = min(abs(alignment.out.Location_time-Spike_times(i)));
                LX(i) = alignment.out.flights(closestIndex(:,i),1);
                LY(i) = alignment.out.flights(closestIndex(:,i),2);
                LZ(i) = alignment.out.flights(closestIndex(:,i),3);
                PH(i) = full(peak_heights(i));
            catch % we need this if Spiketime occurs before/after the location tracking was on..
                continue
            end
        end
        
        % display, in the title, how many bursts there were:
        LX_s = LX;
        LX_s(isnan(LX)) = [];
        disp([num2str(size(LX_s)),' Bursts in flight'])
        PH = mat2gray(PH);
        hold on
        scatter3(LX,LY,LZ,(PH*200)+1,'or','filled');
        title(['Cell no ',num2str(ii),'- ',num2str(size(LX_s)),' Bursts in flight: ' batName ' ' dateSesh ' ' sessionType]);
        xlabel('mm'); ylabel('mm');
    catch % if cell was not active...
        disp('cell not active');
        continue
    end
    
    % Save 'place cells' as jpg and fig files..
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    saveas(gcf,[batName '_' dateSesh '_' sessionType '_placeCell_' num2str(ii) '.tif']);
    savefig(gcf,[batName '_' dateSesh '_' sessionType '_placeCell_' num2str(ii) '.fig']);
    saveas(gcf,[batName '_' dateSesh '_' sessionType '_placeCell_' num2str(ii) '.svg']);

    
    
    
    % Clear the buffer for the next cell:
    clear LX LY LZ closestIndex Spike_times
    
    clf
end

