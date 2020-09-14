pre_dur = 3;     %play with 3-5                                                           %duration of the pre flight period (s):     comparable with flight dur
post_dur = 3;     %play with 3-5                                                           %duration of the post flight period (s):    but shorter than half interflight
dirTop = dir('Ga*');
nDays = [3 7 9]; %which days to look at
%dayCounter = 1; %start the
nRois = [3 11 12];%size(ROIs_gal(day_i,:),2); %number of ROIs
nClusts = 2; %which clusters to look at
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

fig1 = figure();   set(gcf, 'units','normalized','outerposition',[0.2 0 0.5 1]);
sgtitle(['Gal Event Timing Across 3 Days for 3 ROIs: Cluster: 2']);
ha = tight_subplot(3,3,[.06 .03],[.02 .1],[.05 .02]);
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
    
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    % Check if folder exists
    if exist([saveDir1 datestr(now,'yymmdd') filesep batName{dayCounter} '_' dateSesh{dayCounter} '_preDurPostCells'])>0;
        disp('Youve been working today..');
    else
        mkdir([saveDir1 datestr(now,'yymmdd') filesep batName{dayCounter} '_' dateSesh{dayCounter} '_preDurPostCells'])
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep batName{dayCounter} '_' dateSesh{dayCounter} '_preDurPostCells' '\'];
    
   sData = snakeTrace.s; %select data from the S,c_raw, or C matrix
   nFlights = size(flightPaths.clusterIndex{nClusts},1); %find number of flights for that cluster 
   %roiCounter = 1;r
    for roi_i = 1:length(nRois)
        act_pre{day_i,roi_i} = zeros(nFlights,length(sData.normTraceRawPre{nClusts}(1,:,1)));
        act_dur{day_i,roi_i} = zeros(nFlights,length(sData.normTraceRawFlight{nClusts}(1,:,1)));
        act_post{day_i,roi_i} = zeros(nFlights,length(sData.normTraceRawPost{nClusts}(1,:,1)));
        
        for flight_i = 1:nFlights
            act_pre{day_i,roi_i}(flight_i,:) = sData.normTraceRawPre{nClusts}(flight_i,:,ROIs_gal(nDays(day_i),nRois(roi_i)));
            act_dur{day_i,roi_i}(flight_i,:) = sData.normTraceRawFlight{nClusts}(flight_i,:,ROIs_gal(nDays(day_i),nRois(roi_i)));
            act_post{day_i,roi_i}(flight_i,:) = sData.normTraceRawPost{nClusts}(flight_i,:,ROIs_gal(nDays(day_i),nRois(roi_i)));
        end        
        %roiCounter = roiCounter + 1;
        
%         axes(ha(roi_i*day_i*3-2));   imagesc(act_pre{day_i,roi_i});                %set(gca,'xtick',[]);
%         title('Pre');
%         if roi_i ~= 1
%             set(gca,'YTickLabel',[],'XtickLabel',[]);
%         end
        axList = [1 2 3; 4 5 6; 7 8 9];
axes(ha(axList(day_i,roi_i)));   imagesc(act_dur{day_i,roi_i});                %set(gca,'xtick',[]);
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
%         axes(ha(roi_i*day_i*3));    imagesc(act_post{day_i,roi_i});                %set(gca,'xtick',[]);
%         title('Post');
%         if roi_i ~= 1
%             set(gca,'YTickLabel',[],'XtickLabel',[]);
%         end
        hold on;
    
    end
    cd(dirTop(1).folder);
    %dayCounter = dayCounter + 1;
end







% %%
% %converting Angelo's variables to mine
% N = [3 11 12];%size(ROIs_gal(day_i,:),2); %number of ROIs
% %N = size(cellData.results.C_raw,1); %number of ROIs
% traceLength = size(cellData.results.C_raw,2); %length of each neural trace
% CNMFe_Fs = cellData.results.Fs; %imaging sampling rate
% dsFactor =4; %convert from 120hz behavior start frame to 30hz imaging start frame
% 
% %Flight Room references (provvisory)
% xR = +2.85; xL = -2.85; yF = 2.50;  yB = -2.50;  zT = 2.20;                 %Flight volume coordinates
% F1 = [2.56; 1.23; 1.72];    F2 = [2.56; -1.04; 1.72];                       %Feeders coordinates
% F3 = [2.56; 1.43; 0.72];    F4 = [2.56; -1.24; 0.72];                       %Feeders coordinates
% edges_d = {xL:(xR-xL)/10:xR yB:(yF-yB)/10:yF};                              %Edges for density histogram
% 
% %% resample the neural data to account for lost frames
% %Extract variables from Alignment
% img_times = alignment.out.video_times';
% img_sampling = diff(img_times);
% img_sampl_interval = mean(img_sampling);
% img_Fs = 1/img_sampl_interval;
% 
% %Calculate real sampling frequency
% CNMF_Fs_real = round(traceLength/(img_times(end)-img_times(1)));
% CNMF_time = downsample(img_times,round(img_Fs/CNMF_Fs_real));
% 
% %Generate evenly spaced times
% t_even = linspace(img_times(1),img_times(end),length(img_times));
% CNMF_time_even = downsample(t_even,round(img_Fs/CNMF_Fs_real));
% %normalize firing rate by size of neuron
% %A = reshape(cellData.results.A,size(cellData.results.A,2),2); A = permute(A, [3 1 2]);   %Spatial footprints: cell# x pixels x pixels
% %fr_area = sum(A,[2 3])./sum(A,'all');%sum(cellData.results.A,[2 3])./sum(cellData.results.A,'all');                               %Fractional area for each Spatial footprint
% FC_raw = cellData.results.C_raw;%.*fr_area;                                %Raw fluorescence (normalized)
% FC = cellData.results.C;%.*fr_area;                                        %Denoised fluorescence (normalized)
% S = cellData.results.S;%.*fr_area;                                         %Deconvolved spike trace (normalized)
% 
% %Resample at even time intervals
% try
%     FC_raw_rs =    interp1(CNMF_time',FC_raw',CNMF_time_even','linear','extrap')';%cellData.results.C_raw',CNMF_time_even','linear','extrap')';
%     S_rs =     interp1(CNMF_time',S',CNMF_time_even','nearest','extrap')';%cellData.results.S',CNMF_time_even','nearest','extrap')';
% catch
%     recDiff = length(CNMF_time) - length(cellData.results.C_raw);
%     FC_raw_rs =    interp1(CNMF_time(1,1:end-recDiff)',full(cellData.results.C_raw)',CNMF_time_even(1,1:end-recDiff)','linear','extrap')';
%     S_rs =     interp1(CNMF_time(1,1:end-recDiff)',full(cellData.results.S)',CNMF_time_even(1,1:end-recDiff)','nearest','extrap')';
% end
% 
% %Ca-activity and inferred spike-rate are constrained to zero for 5s at the beginning and end of the session and inferred spike-rate is normalized (SD units) and smoothed
% Activity = FC_raw_rs;
% Activity(:,[1:CNMF_Fs_real*5,end-CNMF_Fs_real*5-1:end])=0; 		%cut activity at the start and stop of the recording (5s)
% 
% Rate = normalize(S_rs,2,'scale');       Rate = movmean(Rate,CNMF_Fs_real*0.5,2);        %normalize by STD and smooth on 0.5s
% Rate(:,[1:CNMF_Fs_real*5,end-CNMF_Fs_real*5-1:end])=0; 		%cut activity at the start and stop of the recording (5s)
% 
% until_cluster = 2;%min(length(flightPaths.clusterIndex),100);   %change this if you want less clusters
% 
% %Binning in time and space, p values calculation
% for id_cluster_SI = until_cluster %for each cluster
%     id = [];    %id = find(flight_clus.id==id_cluster_SI); %find all flights that belong to that cluster
%     id = flightPaths.clusterIndex{id_cluster_SI};
%     ROIcounter = 1;
%     for roi_i = N
%         if ~isnan(roi_i)
%             %matrices for keeping the data used in calcuating spatial info
%             bnd_act_pre = zeros(n_bins,size(id,1));
%             bnd_act_dur = zeros(n_bins,size(id,1));
%             bnd_act_pst = zeros(n_bins,size(id,1));
%             Act_pre = zeros(pre_dur*CNMFe_Fs,size(id,1));
%             %Act_dur = zeros(n_bins,size(id,1));
%             Act_pst = zeros(post_dur*CNMFe_Fs,size(id,1));
%             for flight_i=1:size(id,1) %for all flights within the cluster, define the following vectors
%                 %grabbing and binning the velocity and activity along
%                 %with other variables
%                 %convert from behavior time to imaging time
%                 [minValueStart,closestIndexStart] = min(abs(alignment.out.video_times-alignment.out.Location_time(flightPaths.flight_starts_idx(id(flight_i)))));
%                 [minValueEnd,closestIndexEnd] = min(abs(alignment.out.video_times-alignment.out.Location_time(flightPaths.flight_ends_idx(id(flight_i)))));
%                 Act_dur = [];%Act_pre = [];   Act_dur = [];   Act_pst = [];   v_trj = [];
%                 Act_pre(:,flight_i) =  Rate(closestIndexStart-pre_dur*CNMFe_Fs:closestIndexStart-1);
%                 if closestIndexEnd-closestIndexStart == 0 %if flight is 0 frames?
%                     Act_dur =  Rate(closestIndexStart:closestIndexEnd+1);
%                 else
%                     Act_dur =  Rate(closestIndexStart:closestIndexEnd);
%                 end
%                 if closestIndexEnd+post_dur*CNMFe_Fs > length(Rate)
%                     Act_pst(:,flight_i) = Rate(closestIndexEnd+1:length(Rate));
%                 else
%                     Act_pst(:,flight_i) =  Rate(closestIndexEnd+1:closestIndexEnd+post_dur*CNMFe_Fs);
%                 end
%                 
%                 %Temporally binned activity for pre-during-post
%                 flight_dur(1,flight_i) = flightPaths.dur(id(flight_i)); %this comes from the flightPaths output
%                 %bnd_act_pre(:,flight_i) = interp1(linspace(1,100,size(Act_pre,1)),Act_pre,linspace(1,100,n_bins),'linear')';
%                 bnd_act_dur(:,flight_i) = interp1(linspace(1,100,size(Act_dur,2)),Act_dur,linspace(1,100,n_bins),'linear')';
%                 %bnd_act_pst(:,flight_i) = interp1(linspace(1,100,size(Act_pst,2)),Act_pst,linspace(1,100,n_bins),'linear')';
%                 sp_bnd_act(:,flight_i) = interp1(linspace(1,flightPaths.length(id(flight_i)),size(Act_dur,2)),Act_dur,linspace(1,flightPaths.length(id(flight_i)),n_space_bins),'linear')';
%                 
%             end
%             
%             axes(ha(ROIcounter*dayCounter*3-2));   imagesc(Act_pre');                set(gca,'xtick',[]);
%             axes(ha(ROIcounter*dayCounter*3-1));   imagesc(sp_bnd_act');                set(gca,'xtick',[]);
%             axes(ha(ROIcounter*dayCounter*3));    imagesc(Act_pst');                set(gca,'xtick',[]);
%             hold on;
%         end
%         ROIcounter = ROIcounter +1;
%     end
% end
% dayCounter = dayCounter +1;
% %end