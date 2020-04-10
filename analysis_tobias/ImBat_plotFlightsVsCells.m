function [flightVsVelocity,smoothVelocity,smoothAvgSpiking] = ImBat_plotFlightsVsCells(cellData,alignment,flightPaths,analysis_Folder,varargin)

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
smoothAvgSpiking = zscore(smooth(mean((full(cellData.results.C_raw(1:topROILocal,:))),1),25));%-min(data(p).image_data.S(1:numCells,:));


    flightVsVelocity = figure('units','normalized','outerposition',[0 0 1 1]);
    %plot zscored spike data for #cells
    a1 = subplot(topROILocal+11,1,11:topROILocal+11);
    hold on
    for i = 1:topROILocal %[9 31 16 27 26 29 41 46 54 61 68 98 86 87 151]%
        try
            plot(alignment.out.video_timesDS,(zscore(smoothdata(cellData.results.C_raw(i,:),'movmedian',3)))+i*2) %may have to tweak the +i*6 at the end
        catch
            diffVidTimes =  length(cellData.results.C_raw(i,:)) - length(alignment.out.video_timesDS)
            %plot(alignment.out.video_timesDS,(zscore(smoothdata(cellData.results.C_raw(i,1:end-diffVidTimes),'movmedian',3)))+i*2) %may have to tweak the +i*6 at the end
            plot(alignment.out.video_timesDS(1:end-1),(zscore(smoothdata(cellData.results.C_raw(i,:),'movmedian',3)))+i*2) %may have to tweak the +i*6 at the end
        end
    end 
    sgtitle(['Velocity vs Cell Activity: ' batName ' ' dateSesh ' ' sessionType])
    ylabel('z-score dff')
    xlim([0 alignment.out.video_timesDS(end)])%xlim([0 alignment.out.video_timesDS(end)])
    
    %plot smoothed z-scored average firing rate of #cells
    a2 = subplot(topROILocal+11,1,6:9);%topROILocal+3:topROILocal+6);
    hold on
    if length(alignment.out.video_timesDS) == length(smoothAvgSpiking)
        plot(alignment.out.video_timesDS,smoothAvgSpiking)
    else
        diffTimes = length(alignment.out.video_timesDS) - length(smoothAvgSpiking);
        plot(alignment.out.video_timesDS(1:end-diffTimes),smoothAvgSpiking)
    end
    ylim([-1 8])
    xlim([0 alignment.out.video_timesDS(end)])%xlim([0 alignment.out.video_timesDS(end)])
    title('Avg Cell Response')
    ylabel('z-score dff')
    
    %plot the smoothed z-scored velocity of bat & reward vector
    [rewardR,rewardLT,rewardUT] = risetime(alignment.out.RewardVector);
    a3 = subplot(topROILocal+11,1,1:4);%topROILocal+8:topROILocal+11);
    hold on
    plot(alignment.out.Location_time(1:end),smoothVelocity,'color','r')
    plot(alignment.out.Location_time(1:end),(alignment.out.RewardVector-min(alignment.out.RewardVector))/2,'color','b')
    %for ii = 1:length(rewardLT)
       %plot(rewardLT(ii)/120,max(smoothVelocity),'ob');
       %hold on
       %plot((rewardLT(ii)+rewardUT(ii))+10/240,max(smoothVelocity),'or');
       %hold on
       %plot(rewardUT(ii)+5/120,max(smoothVelocity),'ok');
    %end
    title('Velocity')
    ylim([-1 8])
    xlim([0 alignment.out.video_timesDS(end)])
    %xlabel('Time (s)')
    ylabel('z-score velocity')
    
    
    linkaxes([a1, a2, a3], 'x');
    if saveFlag == 1
        saveas(flightVsVelocity, [pwd '\' analysis_Folder '\flights\' label '_flightVsVelocity_handPicked.svg']);
        saveas(flightVsVelocity, [pwd '\' analysis_Folder '\flights\' label '_flightVsVelocity_handPicked.tif']);
        savefig(flightVsVelocity, [pwd '\' analysis_Folder '\flights\' label '_flightVsVelocity_handPicked.fig']);
    end



