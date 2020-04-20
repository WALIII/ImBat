function [plotFiringTrajectory] = ImBat_PlaceCells_Tobias(flightPaths, cellData, alignment,varargin)

offset = 0; % account for slow calcium estimation ~move locations back 100ms in time... This is the knob to turn for 'prospective' coding...
spikeThreshMult = 5; %number of times to multiply std of s vector for determining spike threshold
batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze

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
        case 'loadflag'
            loadFlag = varargin{i+1};
        case 'analysisfolder'
            analysis_Folder = varargin{i+1};
    end
end

%labels for loading and saving data if running independent fromImBat_analyze
if loadFlag == 1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
    %label = [dateSesh '_' sessionType];
    cellData = load([pwd '/processed/Motion_corrected_Data_DS_results.mat']);
    alignment = load([pwd '/processed/Alignment.mat']);
    load([pwd '/' analysis_Folder '/' label '_flightPaths.mat']);
end

% Remove that pesky first flight that starts at 0:0;
try
a = find(alignment.out.flights(:,1) == 0);
a2 = isnan(alignment.out.flights(a(1):end,1));
a3 = find(a2>0);
alignment.out.flights(a(1):a(1)+a3(1),:) = NaN;
alignment.out.Location2 = alignment.out.flights;
catch
end

% Plot the location in space that each cell is active in 1 figure
plotFiringTrajectory =  figure('units','normalized','outerposition',[0 0 0.5 1]);
sgtitle([batName '_' dateSesh '_' sessionType ': Firing Fields']);
for ii = 1:length(cellData.results.S(:,1)); % for each cell
    subplot(ceil(length(cellData.results.S(:,1))/4),4,ii)
    hold on;
    plot(alignment.out.flights(:,1),alignment.out.flights(:,2),'k');% plot the flight trajectory in space
    %plot3(alignment.out.flights(:,1),alignment.out.flights(:,2),alignment.out.flights(:,3),'k');%,'LineWidth',2);% plot the flight trajectory in space

    spikeThresh(ii) = median(cellData.results.S(ii,:)) + std(cellData.results.S(ii,:))*spikeThreshMult; %threshold for eliminating noise from the S matrix
    [~,xy] = find(cellData.results.S(ii,:)>spikeThresh(ii));  % get time neuron is active
    xy = xy(xy<length(alignment.out.video_timesDS));
    Spike_times = alignment.out.video_timesDS(xy)-offset; % convert this to 'spike time'
    peak_heights = cellData.results.S(ii,xy);
    
    LX = zeros(1,length(Spike_times(:,1)));
    LY = zeros(1,length(Spike_times(:,1)));
    LZ = zeros(1,length(Spike_times(:,1)));
    PH = zeros(1,length(Spike_times(:,1)));

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
                disp('cell catch');
                continue
            end
        end
        
        % display, in the title, how many bursts there were:
        LX_s = LX;
        LX_s(isnan(LX)) = [];
        disp([num2str(size(LX_s)),' Bursts in flight'])
        PH = mat2gray(PH);
        hold on
        %scatter(LX,LY,(PH*400)+1,'or','filled');
        scatter(LX,LY,15,'or','filled');
        %scatter3(LX,LY,LZ,(PH*400)+1,'or','filled');
        %uistack(dots,'top');
        % % modify labels for tick marks
        xlim([-3000 3000]);
        ylim([-3000 3000]);        
        xticks = get(gca,'xtick');
        yticks = get(gca,'ytick');
        newlabelsX = arrayfun(@(ax) sprintf('%g', ax/1000), xticks, 'un', 0);
        newlabelsY = arrayfun(@(ay) sprintf('%g', ay/1000), yticks, 'un', 0);
        set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
        xlabel('m'); ylabel('m');
        title(['ROI ' num2str(ii) ' : ' num2str(size(LX_s)) ' Bursts']);  
        hold off
    catch % if cell was not active...
        disp('cell not active');
        continue
    end
    


    
    
    
    
    % Clear the buffer for the next cell:
    clear LX LY LZ closestIndex Spike_times
    
    %clf
end

    % Save 'place cells' as jpg and fig files..
    %set(findall(gcf,'-property','FontSize'),'FontSize',20);
    saveas(gcf,[pwd '\' batName '_' dateSesh '_' sessionType '_placeCell_all.tif']);
    savefig(gcf,[pwd '\placeCells\' batName '_' dateSesh '_' sessionType '_placeCell_all.fig']);
    saveas(gcf,[pwd '\' batName '_' dateSesh '_' sessionType '_placeCell_all.svg']);

    
% Plot the location in space that each cell is active in individual figures
plotFiringTrajectoryIndiv =  figure();
%sgtitle([batName '_' dateSesh '_' sessionType ': Firing Fields']);
for ii = 1:length(cellData.results.S(:,1)); % for each cell
    plot(alignment.out.flights(:,1),alignment.out.flights(:,2),'k');% plot the flight trajectory in space
    %plot3(alignment.out.flights(:,1),alignment.out.flights(:,2),alignment.out.flights(:,3),'k');%,'LineWidth',2);% plot the flight trajectory in space

    spikeThresh(ii) = median(cellData.results.S(ii,:)) + std(cellData.results.S(ii,:))*spikeThreshMult; %threshold for eliminating noise from the S matrix
    [~,xy] = find(cellData.results.S(ii,:)>spikeThresh(ii));  % get time neuron is active
    xy = xy(xy<length(alignment.out.video_timesDS));
    Spike_times = alignment.out.video_timesDS(xy)-offset; % convert this to 'spike time'
    peak_heights = cellData.results.S(ii,xy);
    
    LX = zeros(1,length(Spike_times(:,1)));
    LY = zeros(1,length(Spike_times(:,1)));
    LZ = zeros(1,length(Spike_times(:,1)));
    PH = zeros(1,length(Spike_times(:,1)));

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
                disp('cell catch');
                continue
            end
        end
        
        % display, in the title, how many bursts there were:
        LX_s = LX;
        LX_s(isnan(LX)) = [];
        disp([num2str(size(LX_s)),' Bursts in flight'])
        PH = mat2gray(PH);
        hold on
        %scatter(LX,LY,(PH*400)+1,'or','filled');
        scatter(LX,LY,35,'or','filled');
        %scatter3(LX,LY,LZ,(PH*400)+1,'or','filled');
        %uistack(dots,'top');
        % % modify labels for tick marks
        xlim([-3000 3000]);
        ylim([-3000 3000]);        
        xticks = get(gca,'xtick');
        yticks = get(gca,'ytick');
        newlabelsX = arrayfun(@(ax) sprintf('%g', ax/1000), xticks, 'un', 0);
        newlabelsY = arrayfun(@(ay) sprintf('%g', ay/1000), yticks, 'un', 0);
        set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
        xlabel('m'); ylabel('m');
        title([batName '_' dateSesh '_' sessionType ':ROI ' num2str(ii) ' : ' num2str(size(LX_s)) ' Bursts']);  
        hold off
    catch % if cell was not active...
        disp('cell not active');
        continue
    end
    % Save 'place cells' as jpg and fig files..
    %set(findall(gcf,'-property','FontSize'),'FontSize',20);
    saveas(gcf,[pwd '\' batName '_' dateSesh '_' sessionType '_placeCell_' num2str(ii) '.tif']);
    savefig(gcf,[pwd '\placeCells\' batName '_' dateSesh '_' sessionType '_placeCell_' num2str(ii) '.fig']);
    saveas(gcf,[pwd '\' batName '_' dateSesh '_' sessionType '_placeCell_' num2str(ii) '.svg']);

    % Clear the buffer for the next cell:
    clear LX LY LZ closestIndex Spike_times
    
    clf
end