function [flightVsVelocity,smoothVelocity,smoothAvgSpiking] = ImBat_plotFlightsVsCells(cellData,alignment,flightPaths)

global batName dateSesh sessionType topROI

topROILocal = round((topROI * 0.01 * length(cellData.results.S(:,2)))/2); %top # of cells you want to look at divided by 2 for easier viewing 
%smooth and zscore the velocity
smoothVelocity = zscore(smooth(flightPaths.batSpeed,100));

%image_preprocessing, smooth and average #cells spike responses
smoothAvgSpiking = zscore(smooth(mean((full(cellData.results.S(1:topROILocal,:))),1),50));%-min(data(p).image_data.S(1:numCells,:));

flightVsVelocity = figure('units','normalized','outerposition',[0 0 1 1]);
%plot zscored spike data for #cells
a1 = subplot(topROILocal+11,1,1:topROILocal);
hold on
for i= 1:topROILocal
    plot(alignment.out.video_times,(zscore(cellData.results.S(i,:))/4)+i*6) %may have to tweak the +i*6 at the end
end
title(['Velocity vs Cell Activity: ' batName ' ' dateSesh ' ' sessionType])
ylabel('z-score dff')
xlim([0 alignment.out.video_times(end)])
%plot smoothed z-scored average firing rate of #cells
a2 = subplot(topROILocal+11,1,topROILocal+3:topROILocal+6);
hold on
plot(alignment.out.video_times,smoothAvgSpiking)
ylim([-1 8])
xlim([0 alignment.out.video_times(end)])
title('Avg Cell Response')
ylabel('z-score dff')
%plot the smoothed z-scored velocity of bat
a3 = subplot(topROILocal+11,1,topROILocal+8:topROILocal+11);
hold on
plot(alignment.out.Location_time(1:end),smoothVelocity,'color','r')
title('Velocity')
ylim([-1 8])
xlim([0 alignment.out.video_times(end)])
xlabel('Time (s)');
ylabel('z-score velocity')

    
linkaxes([a1, a2, a3], 'x');