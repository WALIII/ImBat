function [flightVsVelocity,smoothVelocity,smoothAvgSpiking] = ImBat_plotFlightsVsCells(cellData,alignment,flightPaths,varargin)

global topROI

topROI = 60; %top number of ROI to focus on for viewing traces
timeStart = 2620; %starting time of when to focus the plotting limits for zooming in to specific flights/traces
timeEnd = 2680; %ending time of when to focus the plotting limits for zooming in to specific flights/traces
numPlotZoom = 20; %number of traces to view at a time when you zoom in

batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
saveFlag = 0;
%dataSet = 'C_raw';

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
        case 'loadflag'
            loadFlag = varargin{i+1};
    end
end

%labels for loading and saving data if running independent fromImBat_analyze
if loadFlag == 1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
    
    cellData = load([pwd '/processed/Motion_corrected_Data_DS_results.mat']);
    alignment = load([pwd '/processed/Alignment.mat']);
    load([pwd '/analysis/' label '_flightPaths.mat']);
end

if saveFlag ==1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
end

%top # of cells you want to look at divided by 2 for easier viewing
%topROILocal = round((topROI * 0.01 * length(cellData.results.C_raw(:,2)))/2);
%topROILocal = length(cellData.results.C(:,1));
topROILocal = 20;

%smooth and zscore the velocity
smoothVelocity = zscore(smooth(flightPaths.batSpeed,100));

%image_preprocessing, smooth and average #cells spike responses
smoothAvgSpiking = zscore(smooth(mean((full(cellData.results.C_raw(1:topROILocal,:))),1),25));%-min(data(p).image_data.S(1:numCells,:));
%multiply the time limits by imaging sampling rate to go from seconds to frames
timeStartFrames = timeStart*cellData.results.Fs;
timeEndFrames = timeEnd*cellData.results.Fs;
maxPlotNum = floor(length(cellData.results.C(:,1))/numPlotZoom);

for n = 0:maxPlotNum
    flightVsVelocity(n+1) = figure('units','normalized','outerposition',[0 0 0.2 1]);
    %plot zscored spike data for #cells
    a1 = subplot(topROILocal+11,1,11:topROILocal+11);
    hold on
    cellActivityWholeSmooth = zeros(length(cellData.results.C(:,1)),length(cellData.results.C(1,:)));
    for i= (n*numPlotZoom)+1:(n*numPlotZoom)+numPlotZoom%topROILocal
        try
            cellActivityWholeSmooth(i,:) = zscore(smoothdata(cellData.results.C_raw(i,:),'movmedian',3));
            %figure(101);
            %plot(alignment.out.video_times,cellActivityWholeSmooth(i,:)+i*2)
            %hold on
            %figure(102);
            plot(alignment.out.video_times(timeStartFrames-225:timeEndFrames-225),cellActivityWholeSmooth(i,timeStartFrames-225:timeEndFrames-225)+i*2) %may have to tweak the +i*6 at the end
        catch
        end
        %hold on
        %     catch
        %         plot(1:length(cellData.results.C_raw(i,:)),(zscore(smoothdata(cellData.results.C_raw(i,:),'movmedian',8)))+i*2)
        %         disp('Catch plot')
        %     end
    end
    
    title(['Velocity vs Cell Activity: ' batName ' ' dateSesh ' ' sessionType])
    ylabel('z-score dff')
    xlim([timeStart timeEnd])%xlim([0 alignment.out.video_times(end)])
    %plot smoothed z-scored average firing rate of #cells
    a2 = subplot(topROILocal+11,1,6:9);%topROILocal+3:topROILocal+6);
    hold on
    %try
    plot(alignment.out.video_times(timeStartFrames-225:timeEndFrames-225),smoothAvgSpiking(timeStartFrames-225:timeEndFrames-225))
    % catch
    %     plot(1:length(smoothAvgSpiking),smoothAvgSpiking)
    %     disp('catch2')
    % end
    ylim([-1 8])
    xlim([timeStart timeEnd])%xlim([0 alignment.out.video_times(end)])
    title('Avg Cell Response')
    ylabel('z-score dff')
    %plot the smoothed z-scored velocity of bat
    a3 = subplot(topROILocal+11,1,1:4);%topROILocal+8:topROILocal+11);
    hold on
    [minValueStart,closestIndexStart] = min(abs(alignment.out.video_times(timeStartFrames-225)-alignment.out.Location_time));
    [minValueEnd,closestIndexEnd] = min(abs(alignment.out.video_times(timeEndFrames-225)-alignment.out.Location_time));
    plot(alignment.out.Location_time(closestIndexStart:closestIndexEnd),smoothVelocity(closestIndexStart:closestIndexEnd),'color','r')
    %plot(alignment.out.Location_time(1:end),smoothVelocity,'color','r')
    title('Velocity')
    ylim([-1 8])
    xlim([timeStart timeEnd])
    %xlim([0 alignment.out.video_times(end)])
    xlabel('Time (s)')
    ylabel('z-score velocity')
    
    
    linkaxes([a1, a2, a3], 'x');
    if saveFlag == 1
        saveas(flightVsVelocity(n+1), [pwd '\analysis\flights\zoom\' label '_flightVsVelocity' num2str(n+1) '.svg']);
        saveas(flightVsVelocity(n+1), [pwd '\analysis\flights\zoom\' label '_flightVsVelocity' num2str(n+1) '.tif']);
        savefig(flightVsVelocity(n+1), [pwd '\analysis\flights\zoom\' label '_flightVsVelocity' num2str(n+1) '.fig']);
    end
end



