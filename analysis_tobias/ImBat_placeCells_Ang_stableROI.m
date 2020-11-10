function ImBat_placeCells_Ang_stableROI
%function to plot firing fields as red dots against the flight paths of the
%bats for each day focusing only on the stable neurons from ROIs_manual

p_val_analysis = 1; %run the pval analayis?
saveFlag = 1; %do you want to save the figures and output structure?
batId = 'Gen';

if strcmp(batId,'Gal') 
% 15 stable manually selected ROIs across 9 days for Gal
ROIs_manual = [28 20 1 23 12 22 10 8 11 24 NaN 2 21 30 19;
        3 2 10 28 11 1 5 33 8 35 NaN 6 22 32 29;
        4 5 11 24 5 1 16 10 2 18 14 8 25 19 9;
        11 22 4 18 3 1 14 5 19 39 9 17 36 25 8;
        25 6 30 27 3 1 2 9 8 37 NaN 15 31 24 36;
        8 12 35 20 1 10 2 39 9 30 3 14 31 24 11;
        10 7 39 35 3 31 8 22 9 37 5 11 39 17 2;
        9 27 25 45 1 7 8 46 11 33 23 6 42 3 2;
        20 34 29 51 7 10 6 40 16 45 5 8 42 26 43;
        8 2 25 38 16 44 20 7 14 26 3 35 37 24 41;
        1 7 28 NaN 6 17 2 35 16 33 12 11 27 30 34;
        19 4 16 NaN 27 21 3 24 2 29 14 8 26 32 33]; 
g = dir('Ga*');
elseif strcmp(batId,'Gen') 
% 20 stable manually selected ROIs across 5 days for Gen
ROIs_manual = [NaN NaN 10 3 16 12 17 18 27 29 8 9 NaN NaN 21 11 31 15 20 25;
    8 17 5 1 2 6 21 10 18 31 NaN 11 51 53 28 4 38 19 23 20;
    50 54 12 3 48 18 27 15 31 34 NaN NaN 28 NaN 29 25 24 22 38 14;
    9 NaN 7 31 2 22 NaN 20 40 25 13 NaN 34 NaN 26 NaN 45 3 24 21;
    10 36 3 11 2 NaN 21 18 20 9 33 NaN NaN NaN 17 NaN 22 7 30 26]; %14 NaN 3 28 2 6 33 26 18 45 NaN NaN 25 NaN 32 NaN 37 8 28 11
g = dir('Ge*');
end
 % P_value calculation params
    n_bins = 10;     %good around here                                                           %number of bins to divide the pre-during-post flight interval
    n_rep = 1000;    %can lower to 10 for debugging                                                           %number of shufflings
    alfa = 0.05;     %can lower to 0.01 or 0.1                                                           %significance level
    pre_dur = 3;     %play with 3-5                                                           %duration of the pre flight period (s):     comparable with flight dur
    pst_dur = 3;     %play with 3-5                                                           %duration of the post flight period (s):    but shorter than half interflight
    w = gausswin(1); %keep at 1 for no smoothing                                                           %witdh of the gaussian filter (number of bins), use odd values. 1=no filtering
    n_space_bins = 30;  %correlates with space resolution but good around 30 (20cm chunks)                                                        %number of spatial bins
    
%g = dir('Ga*');
z = dir('Z1*');
dirTop = vertcat(g,z); %find all folders in top quality directory

%ROI_duplicate = cell(length(dirTop),1); %make cell for indices of duplicated ROIS


for day_i = 1:length(dirTop)
    %load results data
    try %extract metadata names and enter processed folder
        cd([dirTop(day_i).name filesep 'extracted'])
        flyFolders = dir('*fly*extraction');
        batName{day_i} = flyFolders(end).name(1:3);
        dateSesh{day_i} = flyFolders(end).name(5:10);
        sessionType{day_i} = flyFolders(end).name(12:16);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        if strcmp(batName{day_i}(1),'G')
            cd(dirProcessed(end).name);
        else
            cd(dirProcessed(end).name);
        end
    catch
        cd(dirTop(day_i).name);
        flyFolders = dir('*fly*extraction');
        batName{day_i} = flyFolders(end).name(1:3);
        dateSesh{day_i} = flyFolders(end).name(5:10);
        sessionType{day_i} = flyFolders(end).name(12:16);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        cd(dirProcessed(end).name);
    end
    cellData = load('results.mat');
    alignment = load('Alignment.mat');
    cd(dirProcessed(end).folder);
    
    %load flightPaths data
    %extract metadata names and enter analysis folder
    dirAnalysis = dir('analysis_*');
    if strcmp(batName{day_i}(1),'G')
        cd(dirAnalysis(end).name);
    else
        cd(dirAnalysis(end).name);
    end
    fp = dir('*flightPaths.mat');
    load(fp(end).name); %load flightPaths
    
    close all;
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
% Check if folder exists
if exist([saveDir1 datestr(now,'yymmdd') filesep batName{day_i} '_' dateSesh{day_i} '_preDurPostCells'])>0;
    disp('Youve been working today..');
else
    mkdir([saveDir1 datestr(now,'yymmdd') filesep batName{day_i} '_' dateSesh{day_i} '_preDurPostCells'])
end
saveDir = [saveDir1 datestr(now,'yymmdd') filesep batName{day_i} '_' dateSesh{day_i} '_preDurPostCells' '\'];

    %converting Angelo's variables to mine
    N = size(ROIs_manual(day_i,:),2); %number of ROIs
    %N = size(cellData.results.C_raw,1); %number of ROIs
    T = size(cellData.results.C_raw,2); %length of each neural trace
    CNMFe_Fs = cellData.results.Fs; %imaging sampling rate
    dsFactor =4; %convert from 120hz behavior start frame to 30hz imaging start frame
    
    %Flight Room references (provvisory)
    xR = +2.85; xL = -2.85; yF = 2.50;  yB = -2.50;  zT = 2.20;                 %Flight volume coordinates
    F1 = [2.56; 1.23; 1.72];    F2 = [2.56; -1.04; 1.72];                       %Feeders coordinates
    F3 = [2.56; 1.43; 0.72];    F4 = [2.56; -1.24; 0.72];                       %Feeders coordinates
    edges_d = {xL:(xR-xL)/10:xR yB:(yF-yB)/10:yF};                              %Edges for density histogram
    
    %% resample the neural data to account for lost frames
    %Extract variables from Alignment
    img_times = alignment.out.video_times';
    img_sampling = diff(img_times);
    img_sampl_interval = mean(img_sampling);
    img_Fs = 1/img_sampl_interval;
    
    %Calculate real sampling frequency
    CNMF_Fs_real = round(T/(img_times(end)-img_times(1)));
    CNMF_time = downsample(img_times,round(img_Fs/CNMF_Fs_real));
    
    %Generate evenly spaced times
    t_even = linspace(img_times(1),img_times(end),length(img_times));
    CNMF_time_even = downsample(t_even,round(img_Fs/CNMF_Fs_real));
    %normalize firing rate by size of neuron
    %A = reshape(cellData.results.A,size(cellData.results.A,2),2); A = permute(A, [3 1 2]);   %Spatial footprints: cell# x pixels x pixels
    %fr_area = sum(A,[2 3])./sum(A,'all');%sum(cellData.results.A,[2 3])./sum(cellData.results.A,'all');                               %Fractional area for each Spatial footprint
    FC_raw = cellData.results.C_raw;%.*fr_area;                                %Raw fluorescence (normalized)
    FC = cellData.results.C;%.*fr_area;                                        %Denoised fluorescence (normalized)
    S = cellData.results.S;%.*fr_area;                                         %Deconvolved spike trace (normalized)
    
    %Resample at even time intervals
    try
        FC_raw_rs =    interp1(CNMF_time',FC_raw',CNMF_time_even','linear','extrap')';%cellData.results.C_raw',CNMF_time_even','linear','extrap')';
        S_rs =     interp1(CNMF_time',S',CNMF_time_even','nearest','extrap')';%cellData.results.S',CNMF_time_even','nearest','extrap')';
    catch
        recDiff = length(CNMF_time) - length(cellData.results.C_raw);
        FC_raw_rs =    interp1(CNMF_time(1,1:end-recDiff)',full(cellData.results.C_raw)',CNMF_time_even(1,1:end-recDiff)','linear','extrap')';
        S_rs =     interp1(CNMF_time(1,1:end-recDiff)',full(cellData.results.S)',CNMF_time_even(1,1:end-recDiff)','nearest','extrap')';
    end
    
    %Ca-activity and inferred spike-rate are constrained to zero for 5s at the beginning and end of the session and inferred spike-rate is normalized (SD units) and smoothed
    Activity = FC_raw_rs;
    Activity(:,[1:CNMF_Fs_real*5,end-CNMF_Fs_real*5-1:end])=0; 		%cut activity at the start and stop of the recording (5s)
    
    Rate = normalize(S_rs,2,'scale');       Rate = movmean(Rate,CNMF_Fs_real*0.5,2);        %normalize by STD and smooth on 0.5s
    Rate(:,[1:CNMF_Fs_real*5,end-CNMF_Fs_real*5-1:end])=0; 		%cut activity at the start and stop of the recording (5s)
    
    
    % %Create analysis folder for storing the results
    % analysis_directory=fullfile(pwd,['Analysis_',datestr(now, 'yymmdd_HHMM')]);
    % if ~exist(analysis_directory,'dir')
    %     mkdir(analysis_directory);
    % end
    %% Code in the middle
    
    
    
    %% --------------------------------- REFINED ANALYSIS START FROM HERE------------------------------
    %//////////////////////////////////////////////////////////////////////////////////////////////////
    % Calculations are done on the Rate matrix, which is the S matrix,
    % normalized by standard deviation and smoothed on a 0.5s window
    
    if p_val_analysis
        until_cluster = min(length(flightPaths.clusterIndex),100);   %change this if you want less clusters
        
        %Initialization of matrices and arrays
        frames_to_shift(1)=0;     frames_to_shift(2:n_rep) = randi([10*CNMFe_Fs T-10*CNMFe_Fs],1,n_rep-1);   %Shifting in time (longer than 10s)
        p_val = zeros(3,until_cluster,N); %pre(1),during(2),post(3)/#cluster/#nueron                       %p values for pre, during, post active neurons
        response = zeros(3,until_cluster,N); %sum of each phase                   %Integrated response during pre, during, post periods 'sum(median(Rate)'
        avg_bnd_act = zeros(3*n_bins,until_cluster,N);   %n_bins per section, for each pre/dur/post and for each cluster & cell       %Activity across bins from pre to post
        sp_bnd_response = zeros(n_space_bins,until_cluster,N);  %Spatially binned activity along the trajectory
        sp_bnd_velCel = cell(until_cluster,1);
        S_Info = zeros(2,until_cluster,N); %2 b/c 1dim=actual info, 2dim=p-val                     %bits and p value for spatial information
        %bnd_act_pre = zeros(n_bins,until_cluster,N);
        %bnd_act_post = zeros(n_bins,until_cluster,N);
        
        
        %Binning in time and space, p values calculation
        figure();   set(gcf, 'units','normalized','outerposition',[0.2 0 0.5 1]);
        for id_cluster_SI = 1:3%until_cluster %for each cluster
            id = [];    %id = find(flight_clus.id==id_cluster_SI); %find all flights that belong to that cluster
            id = flightPaths.clusterIndex{id_cluster_SI};
            ROIcounter = 1;
            for cell_n = ROIs_manual(day_i,:)%1:N %for each cell, initialize the below matrices
                if ~isnan(cell_n)
                sgtitle([batName{day_i} ' ' dateSesh{day_i} ' ROI: ' num2str(ROIcounter) ' (' num2str(cell_n) ') Cluster: ' num2str(id_cluster_SI)]);
                disp(['Cell number: ' num2str(ROIcounter) ' (' num2str(cell_n) ')  Trajectory number: ' num2str(id_cluster_SI)]);
                %matrices for keeping the data used in calcuating spatial info
                bnd_act_pre = zeros(n_bins,size(id,1));
                bnd_act_dur = zeros(n_bins,size(id,1));
                bnd_act_pst = zeros(n_bins,size(id,1));
                sp_bnd_act =    zeros(n_space_bins,size(id,1));
                sp_bnd_vel =    zeros(n_space_bins,size(id,1));
                sp_bnd_path =   zeros(n_space_bins,size(id,1));
                sp_bnd_lambda = zeros(n_space_bins,size(id,1));
                spikes = zeros(n_rep,3);
                info = zeros(n_rep,1);
                
                for n = 1:n_rep %for each repetition %repetition 1 is nonshuffled repetition
                    Rate_sh = circshift(Rate(cell_n,:),frames_to_shift(n),2)'; %take rate of cell after doing circular shift by frames_to_shift vector
                    flight_dur = zeros(1,size(id,1));
                    
                    for flight_i=1:size(id,1) %for all flights within the cluster, define the following vectors
                        %grabbing and binning the velocity and activity along
                        %with other variables
                        %convert from behavior time to imaging time
                        [minValueStart,closestIndexStart] = min(abs(alignment.out.video_times-alignment.out.Location_time(flightPaths.flight_starts_idx(id(flight_i)))));
                        [minValueEnd,closestIndexEnd] = min(abs(alignment.out.video_times-alignment.out.Location_time(flightPaths.flight_ends_idx(id(flight_i)))));
                        Act_pre = [];   Act_dur = [];   Act_pst = [];   v_trj = [];
                        Act_pre =  Rate_sh(closestIndexStart-pre_dur*CNMFe_Fs:closestIndexStart-1);
                        if closestIndexEnd-closestIndexStart == 0 %if flight is 0 frames?
                            Act_dur =  Rate_sh(closestIndexStart:closestIndexEnd+1);
                        else
                            Act_dur =  Rate_sh(closestIndexStart:closestIndexEnd);
                        end
                        if closestIndexEnd+pst_dur*CNMFe_Fs > length(Rate_sh)
                            Act_pst = Rate_sh(closestIndexEnd+1:length(Rate_sh));
                        else
                            Act_pst =  Rate_sh(closestIndexEnd+1:closestIndexEnd+pst_dur*CNMFe_Fs);
                        end
                        
                        %Temporally binned activity for pre-during-post
                        flight_dur(1,flight_i) = flightPaths.dur(id(flight_i)); %this comes from the flightPaths output
                        bnd_act_pre(:,flight_i) = interp1(linspace(1,100,size(Act_pre,1)),Act_pre,linspace(1,100,n_bins),'linear')';
                        bnd_act_dur(:,flight_i) = interp1(linspace(1,100,size(Act_dur,1)),Act_dur,linspace(1,100,n_bins),'linear')';
                        bnd_act_pst(:,flight_i) = interp1(linspace(1,100,size(Act_pst,1)),Act_pst,linspace(1,100,n_bins),'linear')';
                        
                        %Spatially binned activity for spatial information
                        %during flight (spatial activity) and interpolating
                        %across the bins
                        v_trj1 = flightPaths.vel(1,~isnan(flightPaths.pos(1,:,id(flight_i))),id(flight_i));
                        v_trj = downsample(v_trj1,4);
                        sp_bnd_act(:,flight_i) = interp1(linspace(1,flightPaths.length(id(flight_i)),size(Act_dur,1)),Act_dur,linspace(1,flightPaths.length(id(flight_i)),n_space_bins),'linear')';
                        sp_bnd_vel(:,flight_i) = interp1(linspace(1,flightPaths.length(id(flight_i)),size(v_trj,2)),v_trj',linspace(1,flightPaths.length(id(flight_i)),n_space_bins),'linear')';
                        sp_bnd_path(:,flight_i) = linspace(1,flightPaths.length(id(flight_i)),n_space_bins)';
                    end
                    
                    %Spatial information (CRITICAL CALCULATION!)
                    prob = 1./mean(sp_bnd_vel,2)*(1/sum(1./mean(sp_bnd_vel,2)));
                    sp_bnd_lambda = sp_bnd_act;
                    lambda =  median(sp_bnd_lambda,2);
                    lambda_ave = lambda'*prob;
                    info(n) = sum(lambda.*prob.*log2((lambda+1e-20)./(lambda_ave+1e-20))); %bits/second
                    
                    %Concatenate activities and plot
                    bnd_act = [bnd_act_pre;bnd_act_dur;bnd_act_pst];
                    subplot(4,4,[1 2 3 5 6 7 9 10 11]);
                    plot(linspace(1,100,3*n_bins),filter(w,1,median(bnd_act,2)),'k');  hold on;
                    if n ==  1 %save specific names for first rep to use for plotting with shaded area
                        avg_bnd_act(:,id_cluster_SI, ROIcounter) = filter(w,1,mean(bnd_act,2));
                        sp_bnd_response(:,id_cluster_SI,ROIcounter) = filter(w,1,lambda);
                        sp_bnd_velCel{id_cluster_SI} = sp_bnd_vel;
                        
                        ciplot(filter(w,1,mean(bnd_act,2))-std(bnd_act,[],2)./sqrt(size(id,1)),filter(w,1,mean(bnd_act,2))+std(bnd_act,[],2)./sqrt(size(id,1)),linspace(1,100,3*n_bins));
                        alpha(0.3);
                        line([33,33], [-max(mean(bnd_act,2)),max(mean(bnd_act,2))],'Color', 'k','LineStyle','--');
                        line([66,66], [-max(mean(bnd_act,2)),max(mean(bnd_act,2))],'Color', 'k','LineStyle','--');
                        xlabel('Flight phase');    ylabel('Activity');
                        
                        subplot(4,4,13);   imagesc(sp_bnd_act');                set(gca,'xtick',[]);
                        subplot(4,4,14);   imagesc(sp_bnd_vel');                set(gca,'xtick',[]);
                        subplot(4,4,15);   imagesc(sp_bnd_vel'.*sp_bnd_act');   set(gca,'xtick',[]);
                    end
                    
                    %Integrated activity in the pre-during-post periods
                    %measure of average firing rate of cell in pre,dur, and
                    %post epochs
                    spikes(n,1) = sum(mean(bnd_act_pre,2))/pre_dur;
                    spikes(n,2) = sum(mean(bnd_act_dur,2))/mean(flight_dur); %uses average flight duration
                    spikes(n,3) = sum(mean(bnd_act_pst,2))/pst_dur;
                end
                
                fig_ord = get(gca,'Children');  set(gca,'Children',circshift(fig_ord,2,1)); hold off;
                %% can quantify all p-values for all spike rates and
                %Pre analysis
                subplot(4,4,4);     hi = histogram(spikes(:,1),'Normalization','pdf'); hold on;
                stem(spikes(1,1),1);   hold off;    title('Pre Spikes');
                p_val(1,id_cluster_SI,ROIcounter) = length(find(spikes(:,1)>=spikes(1,1)))/n_rep; %pval = number of counts that are larger than the spiking of the real trace
                response(1,id_cluster_SI,ROIcounter) = spikes(1,1); %number of responses (spikes) for each flight/average duration of epoch (integrated spike response)
                
                %During analysis
                subplot(4,4,8);     hi = histogram(spikes(:,2),'Normalization','pdf'); hold on;
                stem(spikes(1,2),1);   hold off;    title('During Spikes');
                p_val(2,id_cluster_SI,ROIcounter) = length(find(spikes(:,2)>=spikes(1,2)))/n_rep;
                response(2,id_cluster_SI,ROIcounter) = spikes(1,2);
                
                %Post analysis
                subplot(4,4,12);    hi = histogram(spikes(:,3),'Normalization','pdf'); hold on;
                stem(spikes(1,3),1);   hold off;    title('Post Spikes');
                p_val(3,id_cluster_SI,ROIcounter) = length(find(spikes(:,3)>=spikes(1,3)))/n_rep;
                response(3,id_cluster_SI,ROIcounter) = spikes(1,3);
                
                %Spatial info analysis
                subplot(4,4,16);    hi = histogram(info,'Normalization','pdf'); hold on;
                stem(info(1),1);   hold off;    title('Spatial Info');
                S_Info(1,id_cluster_SI,ROIcounter) = info(1);
                S_Info(2,id_cluster_SI,ROIcounter) = length(find(info>=info(1)))/n_rep;
                
                drawnow();      saveas(gcf,[saveDir batName{day_i} '_' dateSesh{day_i} '_' sessionType{day_i} '_ROI_' num2str(ROIcounter) '_' num2str(cell_n) '_cluster_' num2str(id_cluster_SI) '.tif']);
                end
                ROIcounter = ROIcounter +1;
            end
        end
        
        close all;
    end
    
    %% Look at significant cells
    
    if p_val_analysis
        
        %Putative place cells(**during activity & **spatial info & Peak Firing > 3 average firing)
        %pp_cells is matrix of 0 and 1 to show if cell within each cluster
        %(ncluster x ncells)
        pp_cells = [false(1,N); squeeze(alfa>S_Info(2,2:end,:)>0) & squeeze(alfa>p_val(2,2:end,:)>0) & squeeze(max(sp_bnd_response(:,2:end,:),[],1))>3*mean(Rate,'all')']; %false for cluster #1
        [clusP,roiP] = find(pp_cells==1);
        pp_cells_activity = round(normalize(sp_bnd_response(:,pp_cells),1,'range'),5)'; %normalize activity between 0-1 for putative place fields along 30 space fields
        [sorted_rows,~] = find(pp_cells_activity==1);   sorted_pp_fields = pp_cells_activity(sorted_rows,:);
        figure();   imagesc(sorted_pp_fields,[0 1]);         colormap(viridis); %plots the normalized activity sorted of putative place fields (cells sig/trajectory)
        hold on;
        title([batName{day_i} ' ' dateSesh{day_i} ' ' sessionType{day_i} ' pfields sorted']);
        rname = [];
        for i = 1:length(sorted_rows);
            rname{i} = [num2str(roiP(sorted_rows(i))) '.' num2str(clusP(sorted_rows(i)))];
        end
        set(gca,'ytick',[1:length(sorted_rows)],'yTickLabel',rname);
        ylabel('ROI#.Clust#');
        xlabel('Time Bin');
        saveas(gcf,[saveDir batName{day_i} ' ' dateSesh{day_i} ' ' sessionType{day_i} ' pfields sorted.tif']);
        
        %Putative pre and post cells(**pre/post activity & Peak Firing > 3 average firing)
        %pp_cells is matrix of 0 and 1 to show if cell within each cluster
        %(ncluster x ncells)
        ppre_cells = [false(1,N); squeeze(alfa>p_val(1,2:end,:)>0) & squeeze(max(avg_bnd_act(1:n_bins,2:end,:),[],1))>3*mean(Rate,'all')']; %false for cluster #1
        [clusPre,roiPre] = find(ppre_cells==1);
        ppre_cells_activity = round(normalize(avg_bnd_act(1:n_bins,ppre_cells),1,'range'),5)'; %normalize activity between 0-1 for putative place fields along 30 space fields
        [sorted_rows_pre,~] = find(ppre_cells_activity==1);   sorted_ppre_fields = ppre_cells_activity(sorted_rows_pre,:);
        figure();   imagesc(sorted_ppre_fields,[0 1]);         colormap(viridis); %plots the normalized activity sorted of putative place fields (cells sig/trajectory)
        hold on;
        title([batName{day_i} ' ' dateSesh{day_i} ' ' sessionType{day_i} ' preCells sorted']);
        for i = 1:length(sorted_rows_pre);
            rname{i} = [num2str(roiPre(sorted_rows_pre(i))) '.' num2str(clusPre(sorted_rows_pre(i)))];
        end
        set(gca,'ytick',[1:length(sorted_rows_pre)],'yTickLabel',rname);
        ylabel('ROI#.Clust#');
        xlabel('Time Bin');
        saveas(gcf,[saveDir batName{day_i} ' ' dateSesh{day_i} ' ' sessionType{day_i} ' preCells sorted.tif']);
        
        %post cells
        ppost_cells = [false(1,N); squeeze(alfa>p_val(3,2:end,:)>0) & squeeze(max(avg_bnd_act(n_bins*2+1:n_bins*3,2:end,:),[],1))>3*mean(Rate,'all')']; %false for cluster #1
        [clusPost,roiPost] = find(ppost_cells==1);
        ppost_cells_activity = round(normalize(avg_bnd_act(n_bins*2+1:n_bins*3,ppost_cells),1,'range'),5)'; %normalize activity between 0-1 for putative place fields along 30 space fields
        [sorted_rows_post,~] = find(ppost_cells_activity==1);   sorted_ppost_fields = ppost_cells_activity(sorted_rows_post,:);
        figure();   imagesc(sorted_ppost_fields,[0 1]);         colormap(viridis); %plots the normalized activity sorted of putative place fields (cells sig/trajectory)
        hold on;
        title([batName{day_i} ' ' dateSesh{day_i} ' ' sessionType{day_i} ' postCells sorted']);
        for i = 1:length(sorted_rows_post);
            rname{i} = [num2str(roiPost(sorted_rows_post(i))) '.' num2str(clusPost(sorted_rows_post(i)))];
        end
        set(gca,'ytick',[1:length(sorted_rows_post)],'yTickLabel',rname);
        ylabel('ROI#.Clust#');
        xlabel('Time Bin');
        saveas(gcf,[saveDir batName{day_i} ' ' dateSesh{day_i} ' ' sessionType{day_i} ' postCells sorted.tif']);
       
        %Visualize place fields and calculate centroids
        %outputs every place field as a centroid heatmap along each trajectory
        figure();   set(gcf, 'units','normalized','outerposition',[0.2 0.3 0.5 0.45]);
        [clus,cellNum] = find(pp_cells); %outputs which cells and which clusters have spatial selectivity
        Place_field = struct([]);
        max_position = [];
        for clus_i = 1:length(clus) %for each cluster
            
            %take each flight from a place cell cluster, take the shortest
            %flight, and generate an average trajectory based off the shortest
            %flight
            id = flightPaths.clusterIndex{clus(clus_i)};%find(flight_clus_ds.id==clus(i));
            shortest_flight_i = min(squeeze(sum(~isnan(flightPaths.pos(:,:,id)),2)),[],2);
            ave_trajectory = nanmean(flightPaths.pos(:,1:shortest_flight_i(1),id),3);
            %ave_acceleratn = nanmean(flight_clus_ds.acc(1,1:shortest_flight_i(1),id),3);
            
            npo = size(ave_trajectory,2); %length of average trajectory
            take_off = mean(squeeze(flightPaths.pos(:,1,id)),2);
            map = interp1(linspace(1,100,size(sp_bnd_response(:,clus(clus_i),cellNum(clus_i)),1)),sp_bnd_response(:,clus(clus_i),cellNum(clus_i)),linspace(1,100,npo),'linear')';
            map_color = uint8(round(normalize(map,'range',[2 100]))); %make the mapping between activity and color
            
            %Determine location of the max activity
            [~,max_bin] = max(sp_bnd_velCel{clus(clus_i)}.*sp_bnd_response(:,clus(clus_i),cellNum(clus_i)),[],1,'linear');
            max_position = round(size(ave_trajectory,2)*max_bin/n_space_bins); %find centroid of place field
            %plot the firing activity along the average
            sgtitle([batName{day_i} ' ' dateSesh{day_i} ' ROI: ' num2str(cellNum(clus_i)) ' (' num2str(ROIs_manual(day_i,cellNum(clus_i))) ') Cluster: ' num2str(clus(clus_i))]);
            subplot(121); cmap = viridis(100);
            p = plot3(ave_trajectory(1,:),ave_trajectory(2,:),ave_trajectory(3,:),'r', 'LineWidth',5);
            grid on;
            cdata = [uint8(cmap(map_color,:)*255) uint8(ones(npo,1))].';
            drawnow();  set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cdata);
            hold on;
            textscatter3(take_off(1),take_off(2),take_off(3),"Take-off");
            scatter3(ave_trajectory(1,max_position(1)),ave_trajectory(2,max_position(1)),ave_trajectory(3,max_position(1)));
            xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
            xlabel('x');    ylabel('y');    zlabel('z');
            drawnow();  hold off;
            
            subplot(122);
            plot([1:1:n_space_bins],sp_bnd_response(:,clus(clus_i),cellNum(clus_i)),'LineWidth',3);   xlabel('Space along trajectory');    ylabel('Activity (SD units)');
            xticks([1 n_space_bins]);   xticklabels({'Take-off','Landing'});
            %save info for each place field (position, cell num, clust num)
            Place_field(clus_i).pos = ave_trajectory(:,max_position(1));
            Place_field(clus_i).cell = cellNum(clus_i);
            Place_field(clus_i).clus = clus(clus_i);
            
            saveas(gcf,[saveDir batName{day_i} '_' dateSesh{day_i} '_' sessionType{day_i} '_Place_field' num2str(clus_i) '.tif']);
        end
        close all;
        
        %Calculate percentage of place cells
        perc_place = length(unique(cellNum))./N;
        perc_pre = length(unique(roiPre))./N;
        perc_post = length(unique(roiPost))./N;
        
        %Evaluate distances between centroids for corresponding place fields, within a
        %single cell
        field_dist = [];
        for cell_i = 1:N
            try
            if length(find([Place_field.cell] == cell_i))>1
                cell_pairs = nchoosek(find([Place_field.cell] == cell_i),2);
                for n = 1:size(cell_pairs,1)
                    field_dist = [field_dist; norm([Place_field(cell_pairs(n,1)).pos]-[Place_field(cell_pairs(n,2)).pos])];
                end
            end
            catch
            end
        end
        
        %Calculate percentage of pre, during, post cells on any trajectory (excluding cluster 1)
        percentage = sum(squeeze(any(alfa>p_val(:,2:end,:)>0,2)),2)./N;
        
        %Look at the activity of significantly modulated cells
        %this may have to be corrected because it is not looking at all
        %pre/post/during cells, just looking at during right now
        pt_cells = [false(1,N); squeeze(any(alfa>p_val(:,2:end,:)>0))]; %pt cells have signif pval in any period (pre,during, post)
        [clusAll,roiAll] = find(pt_cells==1);
        pt_cells_activity = round(normalize(avg_bnd_act(:,pp_cells),1,'range'),5)';
        [sorted_rows_all,~] = find(pt_cells_activity==1);   sorted_cells = pt_cells_activity(sorted_rows_all,:);
        figure();       set(gcf, 'units','normalized','outerposition',[0.25 0.25 0.5 0.5]);
        ax1 = subplot(1,3,1);       imagesc(sorted_cells(:,1:n_bins),[0 1]);         colormap(viridis);
        ax2 = subplot(1,3,2);       imagesc(sorted_cells(:,n_bins+1:2*n_bins),[0 1]);
        ax3 = subplot(1,3,3);       imagesc(sorted_cells(:,2*n_bins+1:end),[0 1]);
        ax1.Title.String = 'Pre-Flight';        ax2.Title.String = 'During-Flight'; ax3.Title.String = 'Post-Flight';
        ax1.XLabel.String = 'Bin';              ax2.XLabel.String = 'Bin';          ax3.XLabel.String = 'Bin';
        ax1.YLabel.String = 'ROI#.Clust#';      ax2.YTickLabel= [];                 ax3.YTickLabel= [];
        for i = 1:length(sorted_rows_all);
            rname{i} = [num2str(roiAll(sorted_rows_all(i))) '.' num2str(clusAll(sorted_rows_all(i)))];
        end
        set(gca,'ytick',[1:length(sorted_rows_all)],'yTickLabel',rname);        
        saveas(gcf,[saveDir batName{day_i} '_' dateSesh{day_i} '_' sessionType{day_i} '_activity_sorted.tif']);
        
    end
    
    %% Save (modify this to save only relevant variables)
    if saveFlag == 1
        a_filename = [saveDir batName{day_i} '_' dateSesh{day_i} '_' 'Extracted_trajectories_&_activity',datestr(now, 'yymmdd_HHMM'),'.mat'];
        placeCellsAngStable.meta.n_bins = n_bins;
        placeCellsAngStable.meta.n_rep = n_rep;
        placeCellsAngStable.meta.alfa = alfa;
        placeCellsAngStable.meta.pre_dur = pre_dur;
        placeCellsAngStable.meta.pst_dur = pst_dur;
        placeCellsAngStable.meta.nSpace_bins = n_space_bins;
        placeCellsAngStable.meta.N = N;
        placeCellsAngStable.meta.T = T;
        placeCellsAngStable.meta.CNMFe_Fs = CNMFe_Fs;
        placeCellsAngStable.meta.dsFactor = dsFactor;
        placeCellsAngStable.meta.frames_to_shift = frames_to_shift;
        placeCellsAngStable.meta.batName = batName{day_i};
        placeCellsAngStable.meta.dateSesh = dateSesh{day_i};
        placeCellsAngStable.meta.sessionType = sessionType{day_i};

        placeCellsAngStable.ROIs_manual = ROIs_manual;
        placeCellsAngStable.Rate = Rate;
        placeCellsAngStable.p_val = p_val;
        placeCellsAngStable.response = response;
        placeCellsAngStable.avg_bnd_act = avg_bnd_act;
        placeCellsAngStable.sp_bnd_response = sp_bnd_response;
        placeCellsAngStable.sp_bnd_velCel = sp_bnd_velCel;
        placeCellsAngStable.S_Info = S_Info;
        placeCellsAngStable.max_position = max_position;
        placeCellsAngStable.Place_field = Place_field;
        placeCellsAngStable.percentage = percentage;
        placeCellsAngStable.perc_place = perc_place;
        placeCellsAngStable.perc_pre = perc_pre;
        placeCellsAngStable.perc_post = perc_post;
        placeCellsAngStable.pp_cells = pp_cells;
        placeCellsAngStable.pp_cells_activity = pp_cells_activity;
        placeCellsAngStable.sorted_rows = sorted_rows;
        placeCellsAngStable.sorted_pp_fields = sorted_pp_fields;
        placeCellsAngStable.pt_cells = pt_cells;
        placeCellsAngStable.pt_cells_activity = pt_cells_activity;
        placeCellsAngStable.ppre_cells = ppre_cells;
        placeCellsAngStable.ppre_cells_activity = ppre_cells_activity;
        placeCellsAngStable.sorted_rows_pre = sorted_rows_pre;
        placeCellsAngStable.sorted_ppre_fields = sorted_ppre_fields;
        placeCellsAngStable.ppost_cells = ppost_cells;
        placeCellsAngStable.ppost_cells_activity = ppost_cells_activity;
        placeCellsAngStable.sorted_rows_post = sorted_rows_post;
        placeCellsAngStable.sorted_ppost_fields = sorted_ppost_fields;  
        placeCellsAngStable.clusP = clusP;
        placeCellsAngStable.roiP = roiP;
        placeCellsAngStable.clusPre = clusPre;
        placeCellsAngStable.roiPre = roiPre;
        placeCellsAngStable.clusPost = clusPost;
        placeCellsAngStable.roiPost = roiPost;
        placeCellsAngStable.clusAll = clusAll;
        placeCellsAngStable.roiAll = roiAll;
        placeCellsAngStable.clus = clus;
        placeCellsAngStable.cellNum = cellNum;

        
        save(a_filename,'placeCellsAngStable');
%         if analyze_Ca
%             save([saveDir '/A_', batdate, '.mat'],'A');      %save spatial footprints
%         end
    end
    disp(datestr(now));
    
    cd(dirTop(day_i).folder)
end
close all;