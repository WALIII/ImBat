plotFlag = 0;
saveFlag = 1;
pre_dur = 3;     %play with 3-5                                                           %duration of the pre flight period (s):     comparable with flight dur
post_dur = 3;     %play with 3-5                                                           %duration of the post flight period (s):    but shorter than half interflight
dirTop = dir('Ga*');
nDays = [1:9]; %which days to look at
%dayCounter = 1; %start the
nRois = [1:15];%size(ROIs_gal(day_i,:),2); %number of ROIs
CNMFe_Fs = cellData.results.Fs; %imaging sampling rate
ROIs_gal = [28 20 1 23 12 22 10 8 11 24 NaN 2 21 30 19;
    3 2 10 28 11 1 5 33 8 35 NaN 6 22 32 29;
    4 5 11 24 5 1 16 10 2 18 14 8 25 19 9;
    11 22 4 18 3 1 14 5 19 39 9 17 36 25 8;
    14 3 16 21 2 1 5 7 8 26 NaN 9 27 6 4;
    5 13 41 23 1 21 3 24 6 22 2 25 16 15 7;
    12 3 34 19 2 14 6 15 9 36 5 10 35 20 1;
    25 26 16 32 1 12 4 19 5 28 15 NaN 34 3 2;
    32 34 29 51 7 10 6 40 16 45 5 8 42 26 43];

if plotFlag == 1
    fig1 = figure();   set(gcf, 'units','normalized','outerposition',[0.2 0 0.5 1]);
    sgtitle(['Gal Event Timing Across 3 Days for 3 ROIs: Cluster: 2']);
    ha = tight_subplot(3,3,[.06 .03],[.02 .1],[.05 .02]);
end
%make the matrices to hold the neural data
act_pre = cell(length(nDays),length(nRois));
act_post = cell(length(nDays),length(nRois));
act_dur = cell(length(nDays),length(nRois));

for day_i = 1:length(nDays) %for each day
    %load results data
    try %extract metadata names and enter processed folder
        cd([dirTop(nDays(day_i)).name filesep 'extracted'])
        flyFolders = dir('*fly*extraction');
        batName{dayCounter} = flyFolders(end).name(1:3);
        dateSesh{dayCounter} = flyFolders(end).name(5:10);
        sessionType{dayCounter} = flyFolders(end).name(12:16);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        if strcmp(batName{dayCounter}(1),'G')
            cd(dirProcessed(end).name);
        else
            cd(dirProcessed(end).name); %can change this if need to look at earlier or later processed folders based on batname, date, etc
        end
    catch
        cd(dirTop(nDays(day_i)).name);
        flyFolders = dir('*fly*extraction');
        batName{dayCounter} = flyFolders(end).name(1:3);
        dateSesh{dayCounter} = flyFolders(end).name(5:10);
        sessionType{dayCounter} = flyFolders(end).name(12:16);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        cd(dirProcessed(end).name);
    end
    %cellData = load('results.mat');
    %alignment = load('Alignment.mat');
    cd(dirProcessed(end).folder);
    
    %load flightPaths and snakeTrace data
    %extract metadata names and enter analysis folder
    dirAnalysis = dir('analysis_*');
    if strcmp(batName{dayCounter}(1),'G')
        cd(dirAnalysis(end).name);
    else
        cd(dirAnalysis(end).name);
    end
    fp = dir('*flightPaths.mat');
    load(fp(end).name); %load flightPaths
    sd = dir('*snakePlotData.mat');
    load(sd(end).name);
    
    figh = findall(0,'type','figure');
    other_figures = setdiff(figh, fig1)
    delete(other_figures)
    
    if saveFlag == 1
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    % Check if folder exists
    if exist([saveDir1 datestr(now,'yymmdd') filesep batName{dayCounter} '_' dateSesh{dayCounter} '_preDurPostCells'])>0;
        disp('Youve been working today..');
    else
        mkdir([saveDir1 datestr(now,'yymmdd') filesep batName{dayCounter} '_' dateSesh{dayCounter} '_preDurPostCells'])
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep batName{dayCounter} '_' dateSesh{dayCounter} '_preDurPostCells' '\'];
    end
    
    sData = snakeTrace.cRaw; %select data from the S,c_raw, or C matrix
    nClusts = length(flightPaths.clusterIndex); %which clusters to look at
    nFlights = cell(nClusts,1); 
    lenFlights = cell(nClusts,1);

    
    %for each cluster in nClusts
    for clust_i = 1:nClusts
        nFlights{clust_i} = size(flightPaths.clusterIndex{clust_i},1); %find number of flights for that cluster
        lenFlights{clust_i} = length(sData.normTraceRawPre{clust_i}(1,:,1));
        %roiCounter = 1;r
        for roi_i = 1:length(nRois)
            act_pre{clust_i}{day_i,roi_i} = zeros(nFlights{clust_i},lenFlights{clust_i});
            act_dur{clust_i}{day_i,roi_i} = zeros(nFlights{clust_i},lenFlights{clust_i});
            act_post{clust_i}{day_i,roi_i} = zeros(nFlights{clust_i},lenFlights{clust_i});
            
            for flight_i = 1:nFlights{clust_i}
                if ~isnan(ROIs_gal(nDays(day_i),nRois(roi_i)))
                    act_pre{clust_i}{day_i,roi_i}(flight_i,:) = sData.normTraceRawPre{clust_i}(flight_i,:,ROIs_gal(nDays(day_i),nRois(roi_i)));
                    act_dur{clust_i}{day_i,roi_i}(flight_i,:) = sData.normTraceRawFlight{clust_i}(flight_i,:,ROIs_gal(nDays(day_i),nRois(roi_i)));
                    act_post{clust_i}{day_i,roi_i}(flight_i,:) = sData.normTraceRawPost{clust_i}(flight_i,:,ROIs_gal(nDays(day_i),nRois(roi_i)));
                end
            end
            
            if plotFlag == 1
                axList = [1 2 3; 4 5 6; 7 8 9];
                axes(ha(axList(day_i,roi_i)));   imagesc(act_dur{clust_i}{day_i,roi_i});                %set(gca,'xtick',[]);
                title(['Day ' num2str(nDays(day_i)) ' ROI ' num2str(nRois(roi_i))]);
                xlim([200 300]);
                xt = get(gca,'XTick');
                set(gca,'XTickLabel',round(xt/30,1));
                if roi_i ~= 1
                    set(gca,'YTickLabel',[]);
                else
                    xlabel('Time (s)');
                    ylabel('Flight #');
                end
                hold on;
            end
            
        end
    end
    cd(dirTop(1).folder);
end

activity_allTrials.act_pre = act_pre;
activity_allTrials.act_dur = act_dur;
activity_allTrials.act_post = act_post;
if saveFlag == 1
save(['Gal_200311to200320_activity_allTrials_allClusts' datestr(now,'yymmdd-HHMM') '.mat'],'activity_allTrials');
end

