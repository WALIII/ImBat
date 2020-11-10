batId = 'Gen';
clustNum = 2;
saveFlag = 0;
saveTag = ['sMat'];
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
nDays = size(dataPreDurPost.mean_act_aligned{1},1);
nRois = size(dataPreDurPost.mean_act_aligned{1},2);
lenTrialNeur = size(dataPreDurPost.mean_act_aligned{clustNum}{1,1},2);
lenTrialBehav = size(dataPreDurPost.mean_vel_dur{clustNum}{1,1},2);
%make matrix of correlations between each day
%neural data
%make matrix to hold 3d to get data out of cell
act_allRois = zeros(nDays,nRois,lenTrialNeur);
mean_act_allDay = zeros(nDays,lenTrialNeur);
mean_velXY_dur = zeros(nDays,3,lenTrialBehav);
for day_i = 1:nDays
    for roi_i = 1:nRois
    act_allRois(day_i,roi_i,:) = dataPreDurPost.mean_act_aligned{clustNum}{day_i,roi_i};
    end
    mean_act_allDay(day_i,:) = mean(act_allRois(day_i,:,:));
    %build behavior matrix with velocity, X and Y from just flight times
    mean_velXY_dur(day_i,1,:) = dataPreDurPost.mean_vel_dur{clustNum}{day_i}(1,:); %velocity
    mean_velXY_dur(day_i,2,:) = dataPreDurPost.mean_XY_dur{clustNum}{day_i}(1,1,:); %x position
    mean_velXY_dur(day_i,3,:) = dataPreDurPost.mean_XY_dur{clustNum}{day_i}(1,2,:); %y position  
end

%find correlations between average neural data mean across whole day
for day_i = 1:nDays-1
    for day_ii = 2:nDays
        Rneural = corrcoef(mean_act_allDay(day_i,:),mean_act_allDay(day_ii,:));
        corrNeural(day_i,day_ii) = Rneural(1,2);
        %Rbehav = corrcoef(dataPreDurPost.mean_vel_aligned{clustNum}{day_i}(1,:),datareDurPost.mean_vel_aligned{clustNum}{day_ii}(1,:));
        Rbehav = corrcoef(mean_velXY_dur(day_i,:,:),mean_velXY_dur(day_ii,:,:));
        corrBehav(day_i,day_ii) = Rbehav(1,2);
        %also take correlation for each individual ROI of that day
        for roi_i = 1:nRois 
            RneuralRoi(roi_i,:,:) = corrcoef(act_allRois(day_i,roi_i,:),act_allRois(day_ii,roi_i,:));
            corrNeuralRoi(roi_i,day_i,day_ii) = RneuralRoi(roi_i,1,2);
            RBehavRoi(roi_i,:,:) = corrcoef(mean_velXY_dur(day_i,:,:),mean_velXY_dur(day_ii,:,:));
            corrBehavRoi(roi_i,day_i,day_ii) = RBehavRoi(roi_i,1,2);
        end
    end
end
%plot correlation matrix of neural data from pairwise days
figCorrNeuralBehaviorFootPrint = figure();
sgtitle([batId ' Daily correlations of all mean activity, flights, footprints']);
hold on;
subplot(1,3,1);
imagesc(corrNeural);
title('Neural data');
xlabel('Days');
ylabel('Days');
axis(gca,'equal')
colorbar('southoutside');

%plot correlation matrix of velocity/XY data from pairwise days
subplot(1,3,2);
imagesc(corrBehav);
title('Behavior data');
xlabel('Days');
ylabel('Days');
axis(gca,'equal')
colorbar('southoutside');

if saveFlag == 1
    if strcmp(batId,'Gal')
    saveas(figCorrNeuralBehaviorFootPrint,[saveDir filesep 'Gal_200311to24_corrNeurBehavFootprint_' saveTag '_clust' num2str(clustNum) '_' datestr(now,'yymmdd-HHMM') '.tif']);
    savefig(figCorrNeuralBehaviorFootPrint,[saveDir filesep 'Gal_200311to24_corrNeurBehavFootprint_' saveTag '_clust' num2str(clustNum) '_' datestr(now,'yymmdd-HHMM') '.fig']);
    elseif strcmp(batId,'Gen')
    saveas(figCorrNeuralBehaviorFootPrint,[saveDir filesep 'Gen_200319to24_corrNeurBehavFootprint_' saveTag '_clust' num2str(clustNum) '_' datestr(now,'yymmdd-HHMM') '.tif']);
    savefig(figCorrNeuralBehaviorFootPrint,[saveDir filesep 'Gen_200319to24_corrNeurBehavFootprint_' saveTag '_clust' num2str(clustNum) '_' datestr(now,'yymmdd-HHMM') '.fig']);
       
    end
end



%% make scatter plot of the correlations of behavior and neural data for each ROI
figScatterCorrNeurBehav = figure();
scatter(corrBehavRoi(:,:,:),corrNeuralRoi(:,:,:));
xlabel('XY+V R-val');
ylabel('Dff R-val');
title([batId ' Daily pairwise correlation values across all days & ROIs']);

if saveFlag == 1
    if strcmp(batId,'Gal')
    saveas(figScatterCorrNeurBehav,[saveDir filesep 'Gal_200311to24_corrNeurBehavFootprint_' saveTag '_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.tif']);
    savefig(figScatterCorrNeurBehav,[saveDir filesep 'Gal_200311to24_corrNeurBehavFootprint_' saveTag '_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.fig']);
    elseif strcmp(batId,'Gen')
    saveas(figScatterCorrNeurBehav,[saveDir filesep 'Gen_200319to24_corrNeurBehavFootprint_' saveTag '_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.tif']);
    savefig(figScatterCorrNeurBehav,[saveDir filesep 'Gen_200319to24_corrNeurBehavFootprint_' saveTag '_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.fig']);
       
    end
end
