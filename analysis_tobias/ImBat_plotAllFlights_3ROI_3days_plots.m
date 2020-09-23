saveFlag = 1;
if saveFlag == 1
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    % Check if folder exists
    if exist([saveDir1 datestr(now,'yymmdd')])>0;
        disp('Youve been working today..');
    else
        mkdir([saveDir1 datestr(now,'yymmdd')])
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') '\'];
end

nDays = size(activity_allTrials.act_dur,1); %number of days for the comparison
nRois = size(activity_allTrials.act_dur,2); %number of ROIs/day
minFlightLength = min(cellfun('size',activity_allTrials.act_dur(:,1),2)); %shortest flight length for each day since they may be slightly different lengths
preDur = 90;
postDur = 210;


mean_act_full = cellfun(@mean,activity_allTrials.act_dur,'UniformOutput',false); %take the mean of all the trials for each roi/day
act_aligned = cellfun(@(c) c(:,1:minFlightLength), activity_allTrials.act_dur,'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
mean_act_aligned = cellfun(@mean,act_aligned,'UniformOutput',false); %take the mean of the aligned traces

act_pre = cellfun(@(c) c(:,1:preDur), activity_allTrials.act_dur,'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
mean_act_pre = cellfun(@mean,act_pre,'UniformOutput',false); %take the mean of the aligned traces
act_dur = cellfun(@(c) c(:,preDur+1:minFlightLength-postDur), activity_allTrials.act_dur,'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
mean_act_dur = cellfun(@mean,act_dur,'UniformOutput',false); %take the mean of the aligned traces
act_post = cellfun(@(c) c(:,minFlightLength-postDur+1:minFlightLength), activity_allTrials.act_dur,'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
mean_act_post = cellfun(@mean,act_post,'UniformOutput',false); %take the mean of the aligned traces

%sort the max of each mean
[~,max_mean_act_aligned] = cellfun(@max,mean_act_aligned);
[BMax_mean_act_aligned,IMax_mean_act_aligned] = sort(max_mean_act_aligned,2);
[~,max_mean_act_pre] = cellfun(@max,mean_act_pre);
[BMax_mean_act_pre,IMax_mean_act_pre] = sort(max_mean_act_pre,2);
[~,max_mean_act_dur] = cellfun(@max,mean_act_dur);
[BMax_mean_act_dur,IMax_mean_act_dur] = sort(max_mean_act_dur,2);
[~,max_mean_act_post] = cellfun(@max,mean_act_post);
[BMax_mean_act_post,IMax_mean_act_post] = sort(max_mean_act_post,2);
%% plot all the data sorted by each day's daily peaks
plotSortedDaily = figure();
set(gcf, 'units','normalized','outerposition',[0 0 0.8 0.6]);
sgtitle('Gal Pre/Dur/Post sorted by daily peaks');
ha = tight_subplot(1,9,[.06 .03],[.06 .1],[.05 .02]);
for day_i = 1:nDays
    actCat = cell2mat(mean_act_aligned(day_i,:)');
    axes(ha(day_i));
    imagesc(actCat(IMax_mean_act_aligned(day_i,:),:));
    %set(gca,'xtick',[]);
    title(['Day ' num2str(day_i)]);
    xt = get(gca,'XTick');
    set(gca,'XTickLabel',round(xt/30,1),'YTick',[1:nRois],'yticklabel',IMax_mean_act_aligned(day_i,:));
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
    end
    
end

if saveFlag == 1
    saveas(plotSortedDaily,[saveDir filesep 'Gal_200311to20_snakePlot_sortedDaily' datestr(now,'yymmdd-HHMMSS') '.tif']);
    savefig(plotSortedDaily,[saveDir filesep 'Gal_200311to20_snakePlot_sortedDaily' datestr(now,'yymmdd-HHMMSS') '.fig']);
    
end

%% plot all the data sorted by the activity on day 1
plotSortedDay1 = figure();
set(gcf, 'units','normalized','outerposition',[0 0 0.8 0.6]);
sgtitle('Gal Pre/Dur/Post sorted by day 1');
ha = tight_subplot(1,9,[.06 .03],[.06 .1],[.05 .02]);
for day_i = 1:nDays
    actCat = cell2mat(mean_act_aligned(day_i,:)');
    axes(ha(day_i));
    imagesc(actCat(IMax_mean_act_aligned(1,:),:));
    %set(gca,'xtick',[]);
    title(['Day ' num2str(day_i)]);
    xt = get(gca,'XTick');
    set(gca,'XTickLabel',round(xt/30,1),'YTick',[1:nRois],'yticklabel',IMax_mean_act_aligned(1,:));
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
    end
end

if saveFlag == 1
    saveas(plotSortedDay1,[saveDir filesep 'Gal_200311to20_snakePlot_sortedDay1_cRaw' datestr(now,'yymmdd-HHMMSS') '.tif']);
    savefig(plotSortedDay1,[saveDir filesep 'Gal_200311to20_snakePlot_sortedDay1_cRaw' datestr(now,'yymmdd-HHMMSS') '.fig']);
end

%% plot pre/dur/post stable ROIs based on daily peak
plotSortedDaily_preDurPost = figure();
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
sgtitle('Gal Pre/Dur/Post sorted by daily peak');
ha = tight_subplot(3,9,[.05 .02],[.06 .1],[.05 .02]);
for day_i = 1:nDays
    actCat = cell2mat(mean_act_aligned(day_i,:)'); %extract from cell to matrix
    actCatPre = cell2mat(mean_act_pre(day_i,:)');
    actCatDur = cell2mat(mean_act_dur(day_i,:)');
    actCatPost = cell2mat(mean_act_post(day_i,:)');
    
    axes(ha(day_i));
    imagesc(actCatPre(IMax_mean_act_pre(day_i,:),:));
    %set(gca,'xtick',[]);
    title(['Day ' num2str(day_i) ': Pre']);
    xt = get(gca,'XTick');
    set(gca,'XTickLabel',round(xt/30,1),'YTick',[1:nRois],'yticklabel',IMax_mean_act_pre(1,:));
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
    end
    
    axes(ha(nDays+day_i));
    imagesc(actCatDur(IMax_mean_act_dur(day_i,:),:));
    %set(gca,'xtick',[]);
    title('Dur');
    xt = get(gca,'XTick');
    set(gca,'XTickLabel',round(xt/30,1),'YTick',[1:nRois],'yticklabel',IMax_mean_act_dur(1,:));
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
    end
    
    axes(ha(2*nDays+day_i));
    imagesc(actCatPost(IMax_mean_act_post(day_i,:),:));
    %set(gca,'xtick',[]);
    title(['Post']);
    xt = get(gca,'XTick');
    set(gca,'XTickLabel',round(xt/30,1),'YTick',[1:nRois],'yticklabel',IMax_mean_act_post(1,:));
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
    end
end

if saveFlag == 1
    saveas(plotSortedDay1_preDurPost,[saveDir filesep 'Gal_200311to20_snakePlot_sortedDaily_cRaw_preDurPost' datestr(now,'yymmdd-HHMMSS') '.tif']);
    savefig(plotSortedDay1_preDurPost,[saveDir filesep 'Gal_200311to20_snakePlot_sortedDaily_cRaw_preDurPost' datestr(now,'yymmdd-HHMMSS') '.fig']);
end

%% plot pre/dur/post stable ROIs based on day 1
   plotSortedDay1_preDurPost = figure();
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
sgtitle('Gal Pre/Dur/Post sorted by day 1');
ha = tight_subplot(3,9,[.06 .03],[.06 .1],[.05 .02]);
for day_i = 1:nDays
    actCat = cell2mat(mean_act_aligned(day_i,:)'); %extract from cell to matrix
    actCatPre = cell2mat(mean_act_pre(day_i,:)');
    actCatDur = cell2mat(mean_act_dur(day_i,:)');
    actCatPost = cell2mat(mean_act_post(day_i,:)');
    
    axes(ha(day_i));
    imagesc(actCatPre(IMax_mean_act_pre(1,:),:));
    %set(gca,'xtick',[]);
    title(['Day ' num2str(day_i) ': Pre']);
    xt = get(gca,'XTick');
    set(gca,'XTickLabel',round(xt/30,1),'YTick',[1:nRois],'yticklabel',IMax_mean_act_pre(1,:));
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
    end
    
    axes(ha(nDays+day_i));
    imagesc(actCatDur(IMax_mean_act_dur(1,:),:));
    %set(gca,'xtick',[]);
    title('Dur');
    xt = get(gca,'XTick');
    set(gca,'XTickLabel',round(xt/30,1),'YTick',[1:nRois],'yticklabel',IMax_mean_act_dur(1,:));
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
    end
    
    axes(ha(2*nDays+day_i));
    imagesc(actCatPost(IMax_mean_act_post(1,:),:));
    %set(gca,'xtick',[]);
    title(['Post']);
    xt = get(gca,'XTick');
    set(gca,'XTickLabel',round(xt/30,1),'YTick',[1:nRois],'yticklabel',IMax_mean_act_post(1,:));
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
    end
end

if saveFlag == 1
    saveas(plotSortedDay1_preDurPost,[saveDir filesep 'Gal_200311to20_snakePlot_sortedDay1_cRaw_preDurPost' datestr(now,'yymmdd-HHMMSS') '.tif']);
    savefig(plotSortedDay1_preDurPost,[saveDir filesep 'Gal_200311to20_snakePlot_sortedDay1_cRaw_preDurPost' datestr(now,'yymmdd-HHMMSS') '.fig']);
end