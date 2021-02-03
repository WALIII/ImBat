function [plotFiringTrajectory] = ImBat_PlaceCells_Tobias(flightPaths, cellData, alignment,varargin)

offset = 0; % account for slow calcium estimation ~move locations back 100ms in time... This is the knob to turn for 'prospective' coding...
spikeThreshMult = 5; %number of times to multiply std of s vector for determining spike threshold
speedThresh = 0.7; %threshold for when to eliminate nonflying moments
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
plotFiringTrajectory =  figure('units','normalized','outerposition',[0 0 0.8 1]);
sgtitle([batName ' ' dateSesh  ': Firing Fields']);
ha = tight_subplot(ceil(length(cellData.results.S(:,1))/5),5,[.02 .01],[.01 .08],[.01 .01]);
for ii = 1:length(cellData.results.S(:,1)); % for each cell
    set(0,'CurrentFigure',plotFiringTrajectory);
    axes(ha(ii));
    %a1 = subplot(ceil(length(cellData.results.S(:,1))/5),5,ii);
    hold on;
    plot(alignment.out.flights(:,1),alignment.out.flights(:,2),'k');% plot the flight trajectory in space
    try
        %threshold for eliminating noise from the S matrix
        spikeThresh(ii) = median(cellData.results.S(ii,:)) + std(cellData.results.S(ii,:))*spikeThreshMult;
        [~,xy] = find(cellData.results.S(ii,:)>spikeThresh(ii));  % get time neuron is active
        xy = xy(xy<length(alignment.out.video_times));
        Spike_times = alignment.out.video_times(xy)-offset; % convert this to 'spike time'
        
        %only take data from when bat is flying
        flightVect = alignment.out.flights; %make a new flight vector to eliminate the subthreshold flight_times
        [yz,~] = find(flightPaths.batSpeed<speedThresh); %find when bat is not flying
        flightVect(yz) = NaN; %set nonflying times to NaN
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
                    LX(i) = flightVect(closestIndex(:,i),1);%alignment.out.flights(closestIndex(:,i),1);
                    LY(i) = flightVect(closestIndex(:,i),2);%alignment.out.flights(closestIndex(:,i),2);
                    LZ(i) = flightVect(closestIndex(:,i),3);%alignment.out.flights(closestIndex(:,i),3);
                    if isnan(LX(i))
                        PH(i) = NaN; %make the peak into a NaN for scaling purposes
                    else
                        PH(i) = full(peak_heights(i));
                    end
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
            %scatter(LX,LY,(PH*75)+1,'or','filled');
            axes(ha(ii));
            scatter(LX,LY,(PH*75)+1,'or','filled');
            %scatter3(LX,LY,LZ,(PH*75)+1,'or','filled');
            %uistack(dots,'top');
            % % modify labels for tick marks
            xlim([-3000 3000]);
            ylim([-3000 3000]);
            set(gca,'xticklabel',[],'yticklabel',[]);
            title(['ROI ' num2str(ii) ': ' num2str(size(LX_s)) ' Bursts']);
            %hold off
        catch % if cell was not active...
            disp('cell not active');
            continue
        end
    catch
        xlim([-3000 3000]);
        ylim([-3000 3000]);
        set(gca,'xticklabel',[],'yticklabel',[]);
        title(['ROI ' num2str(n) ': Not Active']);
    end
    
    
    % Plot the location in space that each cell is active in individual figures
    plotFiringTrajectoryIndiv =  figure();
    try
        plot(alignment.out.flights(:,1),alignment.out.flights(:,2),'k');% plot the flight trajectory in space
        hold on;
        disp([num2str(size(LX_s)),' Bursts in flight'])
        PH = mat2gray(PH);
        %scatter(LX,LY,(PH*75)+1,'or','filled');
        scatter(LX,LY,(PH*75)+1,'or','filled');
        %scatter3(LX,LY,LZ,(PH*75)+1,'or','filled');
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
        title([batName ' ' dateSesh ':ROI ' num2str(ii) ': ' num2str(size(LX_s)) ' Bursts']);
        hold off
        
    catch % if cell was not active...
        disp('cell not active');
        xlim([-3000 3000]);
        ylim([-3000 3000]);
        xticks = get(gca,'xtick');
        yticks = get(gca,'ytick');
        newlabelsX = arrayfun(@(ax) sprintf('%g', ax/1000), xticks, 'un', 0);
        newlabelsY = arrayfun(@(ay) sprintf('%g', ay/1000), yticks, 'un', 0);
        set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
        xlabel('m'); ylabel('m');
        title([batName ' ' dateSesh ' ' sessionType ':ROI ' num2str(ii) ': Not Active']);
        hold off
    end
    % Save 'place cells' as jpg and fig files..
    %set(findall(gcf,'-property','FontSize'),'FontSize',20);
    saveas(plotFiringTrajectoryIndiv,[pwd '\' batName '_' dateSesh '_' sessionType '_placeCell_' num2str(ii) '.jpg']);
    savefig(plotFiringTrajectoryIndiv,[pwd filesep batName '_' dateSesh '_' sessionType '_placeCell_' num2str(ii) '.fig']);
   % saveas(plotFiringTrajectoryIndiv,[pwd '\' batName '_' dateSesh '_' sessionType '_placeCell_' num2str(ii) '.svg']);
    
    % Clear the buffer for the next cell:
    clear LX LY LZ closestIndex Spike_times;
    close gcf;
end

% Save 'place cells' as jpg and fig files..
%set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(plotFiringTrajectory,[pwd '\' batName '_' dateSesh '_' sessionType '_placeCell_all.jpg']);
savefig(plotFiringTrajectory,[pwd filesep batName '_' dateSesh '_' sessionType '_placeCell_all.fig']);
%saveas(plotFiringTrajectory,[pwd '\' batName '_' dateSesh '_' sessionType '_placeCell_all.svg']);

close all;