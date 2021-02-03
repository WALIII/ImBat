function [flightVsVelocity,smoothVelocity,smoothAvgSpiking] = ImBat_plotFlightsVsCells(cellData,alignment,flightPaths,varargin)

global topROI

topROI = 60; %top number of ROI to focus on for viewing traces

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
        case 'analysisfolder'
            analysis_Folder = varargin{i+1};
    end
end

%labels for loading and saving data if running independent fromImBat_analyze
if loadFlag == 1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
    processedFolders = dir('processed*');
            processedNewest = sort({processedFolders(end).name});
            processedNewest = char(processedNewest);
    cellData = load([pwd processedNewest '/Motion_corrected_Data_DS_results.mat']);
    alignment = load([pwd processedNewest '/Alignment.mat']);
    load([pwd '/' analysis_Folder '/' label '_flightPaths.mat']);
end

if saveFlag ==1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
end

%top # of cells you want to look at divided by 2 for easier viewing
%topROILocal = round((topROI * 0.01 * length(cellData.results.C_raw(:,2)))/2);
topROILocal = length(cellData.results.C(:,1));


%smooth and zscore the velocity
smoothVelocity = zscore(smooth(flightPaths.batSpeed,100));

%image_preprocessing, smooth and average #cells spike responses
%smoothAvgSpiking = zscore(smooth(mean((full(cellData.results.C_raw(1:topROILocal,:))),1),25));%-min(data(p).image_data.S(1:numCells,:));


    flightVsVelocity = figure('units','normalized','outerposition',[0 0 1 1]);
    %plot zscored spike data for #cells
    a1 = subplot(topROILocal+11,1,11:topROILocal+11);
    hold on
    normSmoothData = zeros(topROILocal,length(cellData.results.C_raw(1,:)));
    for i = 1:topROILocal %[9 31 16 27 26 29 41 46 54 61 68 98 86 87 151]%
        normSmoothData(i,:) = zscore(smoothdata(cellData.results.C_raw(i,:),'movmedian',3)); %normalize each ROI trace
        try
            plot(alignment.out.video_times,(zscore(smoothdata(cellData.results.C_raw(i,:),'movmedian',3)))+i*6) %may have to tweak the +i*6 at the end
        catch %if the tracking recording was turned on before the imaging recording
            try
                diffVidTimes =  length(cellData.results.C_raw(i,:)) - length(alignment.out.video_times);
                plot(alignment.out.video_times,(zscore(smoothdata(cellData.results.C_raw(i,1:end-diffVidTimes),'movmedian',3)))+i*2) %may have to tweak the +i*6 at the end
            catch %if the imaging data misaligned with tracking data by 1 due to dropped frames
                plot(alignment.out.video_times(1:end-1),(zscore(smoothdata(cellData.results.C_raw(i,:),'movmedian',3)))+i*6) %may have to tweak the +i*6 at the end
            end
        end
    end 
    sgtitle(['Velocity vs Cell Activity: ' batName ' ' dateSesh ' ' sessionType])
    ylabel(['z-score dff: ' num2str(length(cellData.results.C(:,1))) ' ROIs'])
    xlabel('Time (s)')
    set(gca,'yticklabel',[]);
    xlim([0 alignment.out.video_times(end)])%xlim([0 alignment.out.video_times(end)])
    
    meanSmoothAvgSpiking = mean(normSmoothData,1); %take average of the normalized smoothed data
    minSmoothAvgSpiking = min(meanSmoothAvgSpiking);
    smoothAvgSpiking = meanSmoothAvgSpiking - minSmoothAvgSpiking;
    
    %plot smoothed z-scored average firing rate of #cells
    a2 = subplot(topROILocal+11,1,6:9);%topROILocal+3:topROILocal+6);
    hold on
    if length(alignment.out.video_times) == length(smoothAvgSpiking)
        plot(alignment.out.video_times,smoothAvgSpiking)
    else
        diffTimes = length(alignment.out.video_times) - length(smoothAvgSpiking);
        plot(alignment.out.video_times(1:end-diffTimes),smoothAvgSpiking)
    end
    ylim([0 max(smoothAvgSpiking)+0.3])
    xlim([0 alignment.out.video_times(end)])%xlim([0 alignment.out.video_times(end)])
    %title('Avg Cell Response')
    ylabel('Avg dff')
    set(gca,'xticklabel',{[]})
    
    %plot the smoothed z-scored velocity of bat & reward vector
    %[rewardR,rewardLT,rewardUT] = risetime(alignment.out.RewardVector);
    a3 = subplot(topROILocal+11,1,1:4);%topROILocal+8:topROILocal+11);
    hold on
    plot(alignment.out.Location_time(1:end),smoothVelocity,'color','r')
    %plot(alignment.out.Location_time(1:end),(alignment.out.RewardVector-abs(min(alignment.out.RewardVector)))/2,'color','b')
%     for ii = 1:length(rewardLT)
%        plot(((alignment.out.Location_time(1)*120)+rewardLT(ii))/120,max(smoothVelocity),'ob');
%        hold on
%        lh = legend('velocity','reward','Location','east');
%        set(lh,'position',[.86 .925 .03 .03]);
%     end
    %title('Velocity')
    ylim([-1 8])
    xlim([0 alignment.out.video_times(end)])
    %xlabel('Time (s)')
    ylabel('v (m/s)')
    set(gca,'xticklabel',{[]})

    linkaxes([a1, a2, a3], 'x');
    if saveFlag == 1
        saveas(flightVsVelocity, [pwd '\' analysis_Folder '\flights\' label '_flightVsVelocity_handPicked.svg']);
        saveas(flightVsVelocity, [pwd '\' analysis_Folder '\flights\' label '_flightVsVelocity_handPicked.jpg']);
        savefig(flightVsVelocity, [pwd '\' analysis_Folder '\flights\' label '_flightVsVelocity_handPicked.fig']);
    end



