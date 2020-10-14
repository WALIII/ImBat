batId = 'Gal';
saveFlag = 1;
cRaw = 0;
clustNum = 2;
Fs_trace = 30;
Fs_behav = 120;
tag ='zscore';

if cRaw == 1
    saveTag = ['cRaw ' tag];
    smoothTrace = 1;
else
    saveTag = ['sMat ' tag];
    smoothTrace = 10;
end
%make saving directory
if saveFlag == 1
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    %saveDir1 = '/Volumes/Tobias_flig/topQualityData/analysis_done/plots/';
    if ~exist([saveDir1 datestr(now,'yymmdd') filesep 'meanAllROIAllDays'])
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'meanAllROIAllDays']);
    else
        disp('You have been working today...');
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'meanAllROIAllDays' filesep];
end

nRois = length(dataPreDurPost.mean_act_aligned{clustNum}(1,:));
lenTrace = length(dataPreDurPost.mean_act_aligned{clustNum}{1,1});
lenBehav = length(dataPreDurPost.mean_vel_aligned{clustNum}{1,1});
nDays = length(dataPreDurPost.mean_act_aligned{clustNum}(:,1));

roiMat = zeros(nDays,nRois,lenTrace);
meanAllRoisEachDay = zeros(nDays,lenTrace);
smoothMeanAllRoisEachDay = zeros(nDays,lenTrace);
stdAllRoisEachDay = zeros(nDays,lenTrace);
meanAllBehavEachDay = zeros(nDays,lenBehav);
smoothMeanAllBehavEachDay = zeros(nDays,lenBehav);
stdAllBehavEachDay = zeros(nDays,lenBehav);
plotMeanAllRoi = figure();
colDays = jet(nDays);
for day_i = 1:nDays
for roi_i = 1:nRois
   roiMat(day_i,roi_i,:) = dataPreDurPost.mean_act_aligned{clustNum}{day_i,roi_i}; 
end

meanAllRoisEachDay(day_i,:) = mean(roiMat(day_i,:,:));
stdAllRoisEachDay(day_i,:) = std(roiMat(day_i,:,:));
meanAllBehavEachDay(day_i,:) = dataPreDurPost.mean_vel_aligned{clustNum}{day_i};
smoothMeanAllRoisEachDay(day_i,:) = zscore(smooth(meanAllRoisEachDay(day_i,:),smoothTrace));
smoothMeanAllRoisEachDay(day_i,:) = smoothMeanAllRoisEachDay(day_i,:) - min(smoothMeanAllRoisEachDay(day_i,:));

subplot(2,1,1)
p1 = plot(1:length(meanAllBehavEachDay(day_i,:)),meanAllBehavEachDay(day_i,:),'Color',colDays(day_i,:),'LineWidth',3);
hold on;
set(gca,'xtick',[]);
ylabel('Vel (m/s)');
xlim([1 lenBehav]);
legInfo{day_i} = (['Day ' num2str(day_i)]);
legend(legInfo);
title('Velocity');

subplot(2,1,2)
p2 = plot(1:length(smoothMeanAllRoisEachDay(day_i,:)),smoothMeanAllRoisEachDay(day_i,:),'Color',colDays(day_i,:),'LineWidth',3);
hold on;
xlim([1 lenTrace]);
title('PSTH');
end
sgtitle([batId ': Mean PSTH of all ROIs for ' num2str(nDays) ' days (' saveTag ')']);
xdat = get(gca,'xtick');
set(gca,'xticklabel',round(xdat/Fs_trace));
xlabel('Time (s)');
if cRaw == 1
    ylabel('Mean df/f');
else
    ylabel('Spikes');
end


if saveFlag == 1
savefig(plotMeanAllRoi,[saveDir batId '_plot_meanAllRoiAllDay_' saveTag '_' datestr(now,'YYmmDD_hhMM') '.fig']);
saveas(plotMeanAllRoi,[saveDir batId '_plot_meanAllRoiAllDay_' saveTag '_' datestr(now,'YYmmDD_hhMM') '.tif']); 
end
