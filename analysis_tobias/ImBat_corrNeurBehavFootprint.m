saveFlag = 0;
saveTag = ['cRaw_smooth'];
if saveFlag == 1
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    % Check if folder exists
    if exist([saveDir1 datestr(now,'yymmdd') filesep 'histogram_peakCorr_acrossDays'])>0;
        disp('Youve been working today..');
    else
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'histogram_peakCorr_acrossDays'])
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'histogram_peakCorr_acrossDays' '\'];
end

nPhases = size(dataPreDurPost.findSelectiveCells,2);
nDays = size(dataPreDurPost.mean_act_aligned,1);
nRois = size(dataPreDurPost.mean_act_aligned,2);
lenTrialNeur = size(dataPreDurPost.mean_act_aligned{1,1},2);
lenTrialBehav = size(dataPreDurPost.mean_vel_aligned{1,1},2);
%make matrix of correlations between each day
%neural data
%make matrix to hold 3d to get data out of cell
act_allRois = zeros(nDays,nRois,lenTrialNeur);
mean_act_allDay = zeros(nDays,lenTrialNeur);
for day_i = 1:nDays
    for roi_i = 1:nRois
    act_allRois(day_i,roi_i,:) = dataPreDurPost.mean_act_aligned{day_i,roi_i};
    end
    mean_act_allDay(day_i,:) = mean(act_allRois(day_i,:,:));
end

%find correlations between average neural data mean across whole day
for day_i = 1:nDays
    for day_ii = 1:nDays
        Rneural = corrcoef(mean_act_allDay(day_i,:),mean_act_allDay(day_ii,:));
        corrNeural(day_i,day_ii) = Rneural(1,2);
        Rbehav = corrcoef(dataPreDurPost.mean_vel_aligned{day_i}(1,:),dataPreDurPost.mean_vel_aligned{day_ii}(1,:));
        corrBehav(day_i,day_ii) = Rbehav(1,2);
    end
end
figCorrNeuralBehaviorFootPrint = figure();
sgtitle('Daily correlations of all mean activity, flights, footprints');
hold on;
subplot(1,3,1);
imagesc(corrNeural);
title('Neural data');
xlabel('Days');
ylabel('Days');
axis(gca,'equal')
colorbar('southoutside');

subplot(1,3,2);
imagesc(corrBehav);
title('Behavior data');
xlabel('Days');
ylabel('Days');
axis(gca,'equal')
colorbar('southoutside');

if saveFlag == 1
    saveas(figCorrNeuralBehaviorFootPrint,[saveDir filesep 'Gal_200311and20_corrNeurBehavFootprint_' saveTag '_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.tif']);
    savefig(figCorrNeuralBehaviorFootPrint,[saveDir filesep 'Gal_200311and20_corrNeurBehavFootprint_' saveTag '_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.fig']);
end
