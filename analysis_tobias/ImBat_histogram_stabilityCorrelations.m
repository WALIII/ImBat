saveFlag = 1;
saveTag = ['cRaw_sMat_smooth'];
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

%set sizes of variables
CNMFe_Fs = 30; %sampling frequency of imaging to divide the timing into seconds
smoothFactor = 10;
nPhases = size(dataPreDurPost.findSelectiveCells,2);
nDays = size(dataPreDurPost.findSelectiveCells,1);
timeCorr = cell(1,nPhases);
timeDiff = cell(1,nPhases);
%create figure
histogram_peakCorr_acrossDays = figure();
sgtitle(['Gal Pairwise Selectively Active Cells Across Days ' saveTag]);
set(gcf, 'units','normalized','outerposition',[0 0 0.8 0.6]);
%find the overlapping ROIs for each day/phase, then gather all the
%differences between max peaks and correlation coefficients for each pair
%of cells that is active across multiple days
for phase_i = 1:nPhases
    overlapCount = 1;
    timeCorr{phase_i} = zeros(nDays,nDays);
    timeDiff{phase_i} = zeros(nDays,nDays);
    for day_i = 1:nDays
        rois = dataPreDurPost.findSelectiveCells{day_i,phase_i};
        nRois = length(dataPreDurPost.findSelectiveCells{day_i,phase_i});
        for roi_i = 1:nRois
            for day_ii = day_i+1:nDays
                w = gausswin(smoothFactor); %smooth the S matrix by gaussian window of 10
                [~,b] = ismember(rois(roi_i),dataPreDurPost.findSelectiveCells{day_ii,phase_i});
                if b > 0
                    %smooth the smatrix and take the difference in peak times
                    if phase_i == 1
                        spikeSmooth1 = filter(w,1,dataPreDurPost.mean_act_pre{day_i,rois(roi_i)});
                        spikeSmooth2 = filter(w,1,dataPreDurPost.mean_act_pre{day_ii,b});
                        timeDiff{phase_i}(day_i,day_ii) = (dataPreDurPost.meanMax_act_pre(day_ii,b) - dataPreDurPost.meanMax_act_pre(day_i,b))/CNMFe_Fs;
                    elseif phase_i == 2
                        spikeSmooth1 = filter(w,1,dataPreDurPost.mean_act_dur{day_i,rois(roi_i)});
                        spikeSmooth2 = filter(w,1,dataPreDurPost.mean_act_dur{day_ii,b});
                        timeDiff{phase_i}(day_i,day_ii) = (dataPreDurPost.meanMax_act_dur(day_ii,b) - dataPreDurPost.meanMax_act_dur(day_i,b))/CNMFe_Fs;
                    elseif phase_i == 3
                        spikeSmooth1 = filter(w,1,dataPreDurPost.mean_act_post{day_i,rois(roi_i)});
                        spikeSmooth2 = filter(w,1,dataPreDurPost.mean_act_post{day_ii,b});
                        timeDiff{phase_i}(day_i,day_ii) = (dataPreDurPost.meanMax_act_post(day_ii,b) - dataPreDurPost.meanMax_act_post(day_i,b))/CNMFe_Fs;
                    end
                    %find correlation coefficient between the two time series
                    R1 = corrcoef(spikeSmooth1,spikeSmooth2);
                    timeCorr{phase_i}(day_i,day_ii) = R1(1,2);
                    overlapCount = overlapCount + 1;
                end
            end
        end
    end
    
    %plot histograms of the time differences and correlations
    phase = ["Pre","Dur","Post"];
    hold on;
    subplot(3,2,2*phase_i-1);
    histogram(timeDiff{phase_i}(find(abs(timeDiff{phase_i})>0)),round(size(timeDiff{phase_i}(find(abs(timeDiff{phase_i})>0)),1)/2)+1,'FaceColor','r');
    try
    xlim([-1*max(abs((timeDiff{phase_i}(find(abs(timeDiff{phase_i})>0)))))-1 max(abs((timeDiff{phase_i}(find(abs(timeDiff{phase_i})>0)))))+1]);
    catch
    end
    if phase_i == 1
        title(['Time diff b/w mean max peaks: ' char(phase(phase_i))]);
    else
        title(char(phase(phase_i)));
    end
    if phase_i ==3
        xlabel('Time (s)');
        ylabel('# ROIs');
     end
    subplot(3,2,2*phase_i);
    histogram(timeCorr{phase_i}(find(abs(timeCorr{phase_i})>0)),round(size(timeCorr{phase_i}(find(abs(timeCorr{phase_i})>0)),1)/2)+1,'FaceColor','b');
    xlim([-1.2 1.2]);
    if phase_i == 1
        title(['CorrCoef of selectively active pairs: ' char(phase(phase_i))]);
    else
        title(char(phase(phase_i)));
    end
    if phase_i == 3
        xlabel('R');
        ylabel('# ROIs');
    end
end

if saveFlag == 1
    saveas(histogram_peakCorr_acrossDays,[saveDir filesep 'Gal_200311and20_histogram_peakDiff_corr_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
    savefig(histogram_peakCorr_acrossDays,[saveDir filesep 'Gal_200311and20_histogram_peakDiff_corr_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
end