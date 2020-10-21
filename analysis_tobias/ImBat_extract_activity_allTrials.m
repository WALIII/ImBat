function [activity_allTrials,flightAligned_vidData] = ImBat_extract_activity_allTrials(batId)
plotFlag = 0;
saveFlag = 1; %do you want to save the figures and output structure?
cRaw = 0;
saveTag = 'allTrials_sMat_large_1to4';
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

if strcmp(batId,'Gal')
    nDays = [1:4]; %which days to look at
    nRois = [1:15];
elseif strcmp(batId,'Gen')
    nDays = [1:5];
    nRois = [1:20];%size(ROIs_manual(day_i,:),2); %number of ROIs
end
%dayCounter = 1; %start the
if strcmp(batId,'Gal') 
% 15 stable manually selected ROIs across 9 days for Gal
ROIs_manual = [28 20 1 23 12 22 10 8 11 24 NaN 2 21 30 19;
    3 2 10 28 11 1 5 33 8 35 NaN 6 22 32 29;
    4 5 11 24 5 1 16 10 2 18 14 8 25 19 9;
    11 22 4 18 3 1 14 5 19 39 9 17 36 25 8;
    14 3 16 21 2 1 5 7 8 26 NaN 9 27 6 4;
    5 13 41 23 1 21 3 24 6 22 2 25 16 15 7;
    12 3 34 19 2 14 6 15 9 36 5 10 35 20 1;
    25 26 16 32 1 12 4 19 5 28 15 NaN 34 3 2;
    32 34 29 51 7 10 6 40 16 45 5 8 42 26 43]; 
dirTop = dir('Ga*');
elseif strcmp(batId,'Gen') 
% 20 stable manually selected ROIs across 5 days for Gen
ROIs_manual = [NaN NaN 10 3 16 12 17 18 27 29 8 9 NaN NaN 21 11 31 15 20 25;
    8 17 5 1 2 6 21 10 18 31 NaN 11 51 53 28 4 38 19 2 20;
    50 54 12 3 48 18 27 15 31 34 NaN NaN 28 NaN 29 25 24 22 38 14;
    8 NaN 4 28 3 18 10 35 42 25 13 NaN 50 39 46 NaN 49 2 32 26;
    14 NaN 3 28 2 6 33 26 18 45 NaN NaN 25 NaN 32 NaN 37 8 28 11];
dirTop = dir('Ge*');
end
if plotFlag == 1
    fig1 = figure();   set(gcf, 'units','normalized','outerposition',[0.2 0 0.5 1]);
    sgtitle(['Gal Event Timing Across 3 Days for 3 ROIs: Cluster: 2']);
    ha = tight_subplot(3,3,[.06 .03],[.02 .1],[.05 .02]);
end
%make the matrices to hold the neural data
act_pre = cell(1,length(nRois));
act_dur = cell(1,length(nRois));
act_post = cell(1,length(nRois));
for clust_ii = 1:length(act_pre)
act_pre{clust_ii} = cell(length(nDays),length(nRois));
act_post{clust_ii} = cell(length(nDays),length(nRois));
act_dur{clust_ii} = cell(length(nDays),length(nRois));
end

for day_i = 1:length(nDays) %for each day
    %load results data
    try %extract metadata names and enter processed folder
        cd([dirTop(nDays(day_i)).name filesep 'extracted'])
        flyFolders = dir('*fly*extraction');
        batName{day_i} = flyFolders(end).name(1:3);
        dateSesh{day_i} = flyFolders(end).name(5:10);
        sessionType{day_i} = flyFolders(end).name(12:16);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        if strcmp(batName{day_i}(1),'G')
            cd(dirProcessed(end).name);
        else
            cd(dirProcessed(end).name); %can change this if need to look at earlier or later processed folders based on batname, date, etc
        end
    catch
        cd(dirTop(nDays(day_i)).name);
        flyFolders = dir('*fly*extraction');
        batName{day_i} = flyFolders(end).name(1:3);
        dateSesh{day_i} = flyFolders(end).name(5:10);
        sessionType{day_i} = flyFolders(end).name(12:16);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        cd(dirProcessed(end).name);
    end
    load('results.mat'); %load cellData
    vidData = load('Motion_corrected_Data.mat'); %load video frame data
    %make gaussian filter based on dff function for spatial filtering
    gSig = 1;
    gSiz = 4.5*gSig; 
    scaling = 1; %scale up the ymax to make less pixelated
    psf = fspecial('gaussian', round(2*gSiz), gSig);
	ind_nonzero = (psf(:)>=max(psf(:,1)));
	psf = psf-mean(psf(ind_nonzero));
	psf(~ind_nonzero) = 0;   % only use pixels within the center disk
%   % Take median, filter, and max of full movie
    Y_med = median(vidData.Y,3);
    Y_min = min(vidData.Y,3);
    Ydff = vidData.Y - Y_min; %subtract min
    Ydff_tFilt = medfilt3(Ydff); %temporal filtering
    Ydff_filt = imfilter(Ydff_tFilt,psf,'symmetric'); %spatial filtering
    YmaxFull{day_i} = max(Ydff_filt,[],3); %take max
    YmaxFull{day_i} = imresize(YmaxFull{day_i},scaling); %resize to eliminate pixelation
    %alignment = load('Alignment.mat');
    cd(dirProcessed(end).folder);
    
    %load flightPaths and snakeTrace data
    %extract metadata names and enter analysis folder
    dirAnalysis = dir('analysis_*');
    if strcmp(batName{day_i}(1),'G')
        cd(dirAnalysis(end).name);
    else
        cd(dirAnalysis(end).name);
    end
    fp = dir('*flightPaths.mat');
    load(fp(end).name); %load flightPaths
    sd = dir('*snakePlotData.mat');
    load(sd(end).name);
    
    if plotFlag == 1
        figh = findall(0,'type','figure');
        other_figures = setdiff(figh, fig1)
        delete(other_figures)
    else
        close all;
    end
    
    if cRaw ==1
        sData = snakeTrace.cRaw; %select data from the s,cRaw, or c matrix
    else
        sData = snakeTrace.s;
    end
    nClusts = length(flightPaths.clusterIndex); %which clusters to look at

    nFlights = cell(nClusts,1); 
    lenFlight = cell(nClusts,1);
    lenPre = cell(nClusts,1);
    lenPost = cell(nClusts,1);
    for clust_i = 1:length(act_pre)
        for roi_i = 1:length(nRois)
            act_pre{clust_i}{day_i,roi_i} = NaN(1,1500);
            act_dur{clust_i}{day_i,roi_i} = NaN(1,1500);
            act_post{clust_i}{day_i,roi_i} = NaN(1,1500);
            vel_pre{clust_i}{day_i} = NaN(1,4000);
            vel_dur{clust_i}{day_i} = NaN(1,4000);
            vel_post{clust_i}{day_i} = NaN(1,4000);
            XY_dur{clust_i}{day_i} = NaN(1,2,4000);
            framesRaw_pre{clust_i}{day_i} = NaN(1,vidData.Ysiz(1),vidData.Ysiz(2),1500);
            framesRaw_dur{clust_i}{day_i} = NaN(1,vidData.Ysiz(1),vidData.Ysiz(2),1500);
            framesRaw_post{clust_i}{day_i} = NaN(1,vidData.Ysiz(1),vidData.Ysiz(2),1500);
            framesFilt_pre{clust_i}{day_i} = NaN(1,vidData.Ysiz(1),vidData.Ysiz(2),1500);
            framesFilt_dur{clust_i}{day_i} = NaN(1,vidData.Ysiz(1),vidData.Ysiz(2),1500);
            framesFilt_post{clust_i}{day_i} = NaN(1,vidData.Ysiz(1),vidData.Ysiz(2),1500);
        end
    end
    %for each cluster in nClusts
    for clust_i = 1:nClusts
        nFlights{clust_i} = size(flightPaths.clusterIndex{clust_i},1); %find number of flights for that cluster
        %find length of pre, dur, post trials for each cluster
        lenPre{clust_i} = length(sData.normTraceRawPre{clust_i}(1,:,1));
        lenFlight{clust_i} = length(sData.normTraceRawFlight{clust_i}(1,:,1));
        lenPost{clust_i} = length(sData.normTraceRawPost{clust_i}(1,:,1));
        lenBehavPre{clust_i} = length(sData.smoothSpeedRawPre{clust_i}(1,:));
        lenBehavDur{clust_i} = length(sData.smoothSpeedRawFlight{clust_i}(1,:));
        lenBehavPost{clust_i} = length(sData.smoothSpeedRawPost{clust_i}(1,:));
        %initialize velocity matrices
        vel_pre{clust_i}{day_i} = zeros(nFlights{clust_i},lenBehavPre{clust_i});
        vel_dur{clust_i}{day_i} = zeros(nFlights{clust_i},lenBehavDur{clust_i});
        vel_post{clust_i}{day_i} = zeros(nFlights{clust_i},lenBehavPost{clust_i});
        XY_dur{clust_i}{day_i} = zeros(nFlights{clust_i},2,lenBehavDur{clust_i});
        framesRaw_pre{clust_i}{day_i} = zeros(nFlights{clust_i},vidData.Ysiz(1),vidData.Ysiz(2),lenPre{clust_i});
        framesRaw_dur{clust_i}{day_i} = zeros(nFlights{clust_i},vidData.Ysiz(1),vidData.Ysiz(2),lenFlight{clust_i});
        framesRaw_post{clust_i}{day_i} = zeros(nFlights{clust_i},vidData.Ysiz(1),vidData.Ysiz(2),lenPost{clust_i});        
        framesFilt_pre{clust_i}{day_i} = zeros(nFlights{clust_i},vidData.Ysiz(1),vidData.Ysiz(2),lenPre{clust_i});
        framesFilt_dur{clust_i}{day_i} = zeros(nFlights{clust_i},vidData.Ysiz(1),vidData.Ysiz(2),lenFlight{clust_i});
        framesFilt_post{clust_i}{day_i} = zeros(nFlights{clust_i},vidData.Ysiz(1),vidData.Ysiz(2),lenPost{clust_i});        
        %roiCounter = 1;r
        for roi_i = 1:length(nRois)
            %if nFlights{clust_i} > 0
                act_pre{clust_i}{day_i,roi_i} = zeros(nFlights{clust_i},lenPre{clust_i});
                act_dur{clust_i}{day_i,roi_i} = zeros(nFlights{clust_i},lenFlight{clust_i});
                act_post{clust_i}{day_i,roi_i} = zeros(nFlights{clust_i},lenPost{clust_i});
            %else
                
            %end
                
            for flight_i = 1:nFlights{clust_i}
                if ~isnan(ROIs_manual(nDays(day_i),nRois(roi_i)))
                    act_pre{clust_i}{day_i,roi_i}(flight_i,:) = sData.normTraceRawPre{clust_i}(flight_i,:,ROIs_manual(nDays(day_i),nRois(roi_i)));
                    act_dur{clust_i}{day_i,roi_i}(flight_i,:) = sData.normTraceRawFlight{clust_i}(flight_i,:,ROIs_manual(nDays(day_i),nRois(roi_i)));
                    act_post{clust_i}{day_i,roi_i}(flight_i,:) = sData.normTraceRawPost{clust_i}(flight_i,:,ROIs_manual(nDays(day_i),nRois(roi_i)));
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
        Fs_cnmfe = 30;
        prePad = 3*Fs_cnmfe;
        postPad = 7 * Fs_cnmfe;
        for flight_i = 1:nFlights{clust_i}
            %build matrices with the imaging raw frames for flight-aligned max projections
            %cut out raw data, subtract min, temporal and spatial filtering
            framesRaw_pre{clust_i}{day_i}(flight_i,:,:,:) = vidData.Y(:,:,sData.startIdxTrace{clust_i}(flight_i)-sData.preFlightPadCalcium:sData.startIdxTrace{clust_i}(flight_i));
            Y_minPre = min(framesRaw_pre{clust_i}{day_i}(flight_i,:,:,:),4);
            framesDff_pre = framesRaw_pre{clust_i}{day_i}(flight_i,:,:,:) - Y_minPre;
            framesDff_pre = squeeze(framesDff_pre);
            framesTFilt_pre = medfilt3(framesDff_pre);
            framesFilt_pre{clust_i}{day_i}(flight_i,:,:,:) = imfilter(framesTFilt_pre,psf,'symmetric');
            try
            framesRaw_dur{clust_i}{day_i}(flight_i,:,:,:) = vidData.Y(:,:,sData.startIdxTrace{clust_i}(flight_i)-prePad:...
                sData.endIdxTrace{clust_i}(flight_i) + (sData.medDur{clust_i}-sData.dur{clust_i}(flight_i))+ postPad);
            framesRaw_post{clust_i}{day_i}(flight_i,:,:,:) = vidData.Y(:,:,sData.endIdxTrace{clust_i}(flight_i):sData.endIdxTrace{clust_i}(flight_i)+sData.postFlightPadCalcium);
            catch
            framesRaw_dur{clust_i}{day_i}(flight_i,:,:,1:length(vidData.Y(:,:,sData.startIdxTrace{clust_i}(flight_i)-prePad:end)))...
                = vidData.Y(:,:,sData.startIdxTrace{clust_i}(flight_i)-prePad:end);
            framesRaw_post{clust_i}{day_i}(flight_i,:,:,1:length(vidData.Y(:,:,sData.endIdxTrace{clust_i}(flight_i):end)))...
                = vidData.Y(:,:,sData.endIdxTrace{clust_i}(flight_i):end);
            end
            %cut out raw data, subtract min, temporal and spatial filtering
            Y_minDur = min(framesRaw_dur{clust_i}{day_i}(flight_i,:,:,:),3);
            framesDff_dur = framesRaw_dur{clust_i}{day_i}(flight_i,:,:,:) - Y_minDur;
            framesDff_dur = squeeze(framesDff_dur);
            framesTFilt_dur = medfilt3(squeeze(framesDff_dur));
            framesFilt_dur{clust_i}{day_i}(flight_i,:,:,:) = imfilter(framesTFilt_dur,psf,'symmetric');
            Y_minPost = min(framesRaw_post{clust_i}{day_i}(flight_i,:,:,:),3);
            framesDff_post = framesRaw_post{clust_i}{day_i}(flight_i,:,:,:) - Y_minPost;
            framesDff_post = squeeze(framesDff_post);
            framesTFilt_post = medfilt3(squeeze(framesDff_post));
            framesFilt_post{clust_i}{day_i}(flight_i,:,:,:) = imfilter(framesTFilt_post,psf,'symmetric');

            %add data to matrices for velocity of the active trials
            if ~isnan(ROIs_manual(nDays(day_i),nRois(roi_i)))
                vel_pre{clust_i}{day_i}(flight_i,:) = sData.smoothSpeedRawPre{clust_i}(flight_i,:);
                vel_dur{clust_i}{day_i}(flight_i,:) = sData.smoothSpeedRawFlight{clust_i}(flight_i,:);
                vel_post{clust_i}{day_i}(flight_i,:) = sData.smoothSpeedRawPost{clust_i}(flight_i,:);
                XY_dur{clust_i}{day_i}(flight_i,1,:) = sData.posFlight{clust_i}(flight_i,1,:); %flightPaths.pos(1,:,flightPaths.clusterIndex{clust_i}(flight_i));
                XY_dur{clust_i}{day_i}(flight_i,2,:) = sData.posFlight{clust_i}(flight_i,2,:);            %flightPaths.pos(2,:,flightPaths.clusterIndex{clust_i}(flight_i));
            end
        end
        %take average across all trials to get 'average' flight activity
        meanFrames_pre{clust_i}{day_i} = squeeze(nanmean(framesFilt_pre{clust_i}{day_i},1));
        meanFrames_dur{clust_i}{day_i} = squeeze(nanmean(framesFilt_dur{clust_i}{day_i},1));
        meanFrames_post{clust_i}{day_i} = squeeze(nanmean(framesFilt_post{clust_i}{day_i},1));
        %find max projection of this mean
        maxMeanFrames_pre{clust_i}{day_i} = max(meanFrames_pre{clust_i}{day_i},[],3);
        maxMeanFrames_pre{clust_i}{day_i} = imresize(maxMeanFrames_pre{clust_i}{day_i},scaling);
        maxMeanFrames_dur{clust_i}{day_i} = max(meanFrames_dur{clust_i}{day_i},[],3);
        maxMeanFrames_dur{clust_i}{day_i} = imresize(maxMeanFrames_dur{clust_i}{day_i},scaling);
        maxMeanFrames_post{clust_i}{day_i} = max(meanFrames_post{clust_i}{day_i},[],3);        
        maxMeanFrames_post{clust_i}{day_i} = imresize(maxMeanFrames_post{clust_i}{day_i},scaling);
    end
    cd(dirTop(1).folder);
    disp(num2str(day_i))
    clear vidData;
end


activity_allTrials.batId = batId;
activity_allTrials.act_pre = act_pre;
activity_allTrials.act_dur = act_dur;
activity_allTrials.act_post = act_post;
activity_allTrials.vel_pre = vel_pre;
activity_allTrials.vel_dur = vel_dur;
activity_allTrials.vel_post = vel_post;
activity_allTrials.XY_dur = XY_dur;
activity_allTrials.meanFrames_pre = meanFrames_pre;
activity_allTrials.meanFrames_dur = meanFrames_dur;
activity_allTrials.meanFrames_post = meanFrames_post;
activity_allTrials.maxMeanFrames_pre = maxMeanFrames_pre;
activity_allTrials.maxMeanFrames_dur = maxMeanFrames_dur;
activity_allTrials.maxMeanFrames_post = maxMeanFrames_post;
activity_allTrials.YmaxFull = YmaxFull;

%flightAligned_vidData.frames_pre = frames_pre;
flightAligned_vidData.framesRaw_dur = framesRaw_dur;
flightAligned_vidData.framesFilt_dur = framesFilt_dur;
%flightAligned_vidData.frames_post = frames_post;

if saveFlag == 1
    if strcmp(batId,'Gal')
        save([saveDir 'Gal_200311to200320_activity_allTrials_allClusts_' saveTag '.mat'],'activity_allTrials');
        save([saveDir 'Gal_200311to200320_flightAligned_vidData_' saveTag '.mat'],'flightAligned_vidData','-v7.3');        
    elseif strcmp(batId,'Gen')
        save([saveDir 'Gen_200319to200324_activity_allTrials_allClusts_' saveTag '.mat'],'activity_allTrials');    
        save([saveDir 'Gen_200319to200324_flightAligned_vidData_' saveTag '.mat'],'flightAligned_vidData','-v7.3');    
    end
end

