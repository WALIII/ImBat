function ImBat_placeCells_Ang(flightPaths, cellData, alignment,varargin)

batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
p_val_analysis = 1;
save_data = 1;

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'batname'
            batName=varargin{i+1};
        case 'datesesh'
            dateSesh = varargin{i+1};
        case 'sessiontype'
            sessionType = varargin{i+1};
        case 'loadflag'
            loadFlag = varargin{i+1};
        case 'analysisfolder'
            analysis_Folder = varargin{i+1};
    end
end

%labels for loading and saving data if running independent fromImBat_analyze
if loadFlag == 1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
    %label = [dateSesh '_' sessionType];
    cellData = load([pwd '/processed/Motion_corrected_Data_DS_results.mat']);
    alignment = load([pwd '/processed/Alignment.mat']);
    load([pwd '/' analysis_Folder '/' label '_flightPaths.mat']);
end

% P_value calculation params
n_bins = 30;     %good around here                                                           %number of bins to divide the pre-during-post flight interval
n_rep = 500;    %can lower to 10 for debugging                                                           %number of shufflings
alfa = 0.05;     %can lower to 0.01 or 0.1                                                           %significance level
pre_dur = 3;     %play with 3-5                                                           %duration of the pre flight period (s):     comparable with flight dur
pst_dur = 3;     %play with 3-5                                                           %duration of the post flight period (s):    but shorter than half interflight
w = gausswin(1); %keep at 1 for no smoothing                                                           %witdh of the gaussian filter (number of bins), use odd values. 1=no filtering
n_space_bins = 30;  %correlates with space resolution but good around 30 (20cm chunks)                                                        %number of spatial bins

%converting Angelo's variables to mine
N = size(cellData.results.C_raw,1); %number of ROIs
T = size(cellData.results.C_raw,2); %length of each neural trace
CNMFe_Fs = cellData.results.Fs; %imaging sampling rate
%dsFactor =4; %convert from 120hz behavior start frame to 30hz imaging start frame

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
S = full(cellData.results.S);%.*fr_area;                                         %Deconvolved spike trace (normalized)

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
    until_cluster = 4;%min(length(flightPaths.clusterIndex),100);   %change this if you want less clusters
    
    %Create folder to store figures
    figures_directory1=fullfile(pwd,'Spatial_information');
    delete([figures_directory1 '\*']);
    if exist(figures_directory1,'dir')~=7
        mkdir(figures_directory1);
    end
    
    %Initialization of matrices and arrays
    frames_to_shift(1)=0;     frames_to_shift(2:n_rep) = randi([10*CNMFe_Fs T-10*CNMFe_Fs],1,n_rep-1);   %Shifting in time (longer than 10s)
    p_val = zeros(3,until_cluster,N); %pre(1),during(2),post(3)/#cluster/#nueron                       %p values for pre, during, post active neurons
    response = zeros(3,until_cluster,N); %sum of each phase                   %Integrated response during pre, during, post periods 'sum(median(Rate)'
    avg_bnd_act = zeros((3*n_bins)-3,until_cluster,N);   %n_bins per section, for each pre/dur/post and for each cluster & cell       %Activity across bins from pre to post
    sp_bnd_response = zeros(n_space_bins,until_cluster,N);  %Spatially binned activity along the trajectory
    sp_bnd_velCel = cell(until_cluster,1);
    S_Info = zeros(2,until_cluster,N); %2 b/c 1dim=actual info, 2dim=p-val                     %bits and p value for spatial information
    
    %Binning in time and space, p values calculation
    figure();   set(gcf, 'units','normalized','outerposition',[0.2 0 0.5 1]);
    for id_cluster_SI = 2:until_cluster %for each cluster
        id = [];    %id = find(flight_clus.id==id_cluster_SI); %find all flights that belong to that cluster
        id = flightPaths.clusterIndex{id_cluster_SI};
        for cell_n = 1:N %for each cell, initialize the below matrices
            sgtitle(['ROI: ' num2str(cell_n) ' Cluster: ' num2str(id_cluster_SI)]);
            disp(['Cell number: ' num2str(cell_n) '  Trajectory number: ' num2str(id_cluster_SI)]);
            %matrices for keeping the data used in calcuating spatial info
            bnd_act_pre = zeros(n_bins,size(id,1));
            bnd_act_dur = zeros(n_bins,size(id,1));
            bnd_act_pst = zeros(n_bins,size(id,1));
            sp_bnd_act =    zeros(n_space_bins,size(id,1));
            sp_bnd_act_pre = zeros(n_space_bins,size(id,1));
            sp_bnd_act_pst = zeros(n_space_bins,size(id,1));
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
                    [minValueStart,closestIndexStart] = min(abs(CNMF_time_even-alignment.out.Location_time(flightPaths.flight_starts_idx(id(flight_i)))));%min(abs(alignment.out.video_times-alignment.out.Location_time(flightPaths.flight_starts_idx(id(flight_i)))));
                    [minValueEnd,closestIndexEnd] = min(abs(CNMF_time_even-alignment.out.Location_time(flightPaths.flight_ends_idx(id(flight_i)))));%min(abs(alignment.out.video_times-alignment.out.Location_time(flightPaths.flight_ends_idx(id(flight_i)))));
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
                    sp_bnd_act_pre(:,flight_i) = interp1(linspace(1,size(Act_pre,1),size(Act_pre,1)),Act_pre,linspace(1,size(Act_pre,1),n_space_bins),'linear')';
                    sp_bnd_act_pst(:,flight_i) = interp1(linspace(1,size(Act_pst,1),size(Act_pst,1)),Act_pst,linspace(1,size(Act_pst,1),n_space_bins),'linear')';
                    
                    sp_bnd_vel(:,flight_i) = interp1(linspace(1,flightPaths.length(id(flight_i)),size(v_trj,2)),v_trj',linspace(1,flightPaths.length(id(flight_i)),n_space_bins),'linear')';
                    sp_bnd_path(:,flight_i) = linspace(1,flightPaths.length(id(flight_i)),n_space_bins)';
                end
                
                %Spatial information (CRITICAL CALCULATION!)
                prob = 1./mean(sp_bnd_vel,2)*(1/sum(1./mean(sp_bnd_vel,2)));
                sp_bnd_lambda = sp_bnd_act;
                lambda =  median(sp_bnd_lambda,2);
                lambda_pre = median(sp_bnd_act_pre,2);
                lambda_pst = median(sp_bnd_act_pst,2);
                lambda_ave = lambda'*prob;
                info(n) = sum(lambda.*prob.*log2((lambda+1e-20)./(lambda_ave+1e-20))); %bits/second
                
                %Concatenate activities and plot
                bnd_act = [bnd_act_pre(1:end-1,:);bnd_act_dur(1:end-1,:);bnd_act_pst(1:end-1,:)];
                subplot(4,4,[1 2 3 5 6 7 9 10 11]);
                plot(linspace(1,90,(3*n_bins)-3),filter(w,1,median(bnd_act,2)),'k');  hold on;
                if n ==  1 %save specific names for first rep to use for plotting with shaded area
                    avg_bnd_act(:,id_cluster_SI, cell_n) = filter(w,1,median(bnd_act,2));
                    sp_bnd_response(:,id_cluster_SI,cell_n) = filter(w,1,lambda);
                    sp_bnd_response_pre(:,id_cluster_SI,cell_n) = filter(w,1,lambda_pre);
                    sp_bnd_response_pst(:,id_cluster_SI,cell_n) = filter(w,1,lambda_pst);
                    
                    sp_bnd_velCel{id_cluster_SI} = sp_bnd_vel;
                    %standard error mean for bounded line
                    ciplot(filter(w,1,median(bnd_act,2))-std(bnd_act,[],2)./sqrt(size(id,1)),filter(w,1,median(bnd_act,2))+std(bnd_act,[],2)./sqrt(size(id,1)),linspace(1,90,(3*n_bins)-3));
                    alpha(0.3);
                    line([30,30], [-max(median(bnd_act,2)),max(median(bnd_act,2))],'Color', 'k','LineStyle','--');
                    line([60,60], [-max(median(bnd_act,2)),max(median(bnd_act,2))],'Color', 'k','LineStyle','--');
                    xlabel('Flight phase');    ylabel('Activity');
                    
                    subplot(4,4,13);   imagesc([sp_bnd_act_pre; sp_bnd_act; sp_bnd_act_pst]');  title('PreDurPost Raster');              %set(gca,'xtick',[]);
                    subplot(4,4,14);   imagesc(sp_bnd_vel');  title('Velocity');               %set(gca,'xtick',[]);
                    subplot(4,4,15);   imagesc(sp_bnd_vel'.*sp_bnd_act'); title('Normalized Flight Raster');   %set(gca,'xtick',[]);
                end
                
                %Integrated activity in the pre-during-post periods
                %measure of average firing rate of cell in pre,dur, and
                %post epochs
                spikes(n,1) = sum(median(bnd_act_pre,2))/pre_dur; %mean vs median
                spikes(n,2) = sum(median(bnd_act_dur,2))/mean(flight_dur); %uses average flight duration
                spikes(n,3) = sum(median(bnd_act_pst,2))/pst_dur;
            end
            
            fig_ord = get(gca,'Children');  set(gca,'Children',circshift(fig_ord,2,1)); hold off;
            %% can quantify all p-values for all spike rates and 
            %Pre analysis
            subplot(4,4,4);     hi = histogram(spikes(:,1),'Normalization','pdf'); hold on;
            stem(spikes(1,1),1);   hold off;    title('Pre Spikes');
            p_val(1,id_cluster_SI,cell_n) = length(find(spikes(:,1)>=spikes(1,1)))/n_rep; %pval = number of counts that are larger than the spiking of the real trace
            response(1,id_cluster_SI,cell_n) = spikes(1,1); %number of responses (spikes) for each flight/average duration of epoch (integrated spike response)
            
            %During analysis
            subplot(4,4,8);     hi = histogram(spikes(:,2),'Normalization','pdf'); hold on;
            stem(spikes(1,2),1);   hold off;    title('During Spikes');
            p_val(2,id_cluster_SI,cell_n) = length(find(spikes(:,2)>=spikes(1,2)))/n_rep;
            response(2,id_cluster_SI,cell_n) = spikes(1,2);
            
            %Post analysis
            subplot(4,4,12);    hi = histogram(spikes(:,3),'Normalization','pdf'); hold on;
            stem(spikes(1,3),1);   hold off;    title('Post Spikes');
            p_val(3,id_cluster_SI,cell_n) = length(find(spikes(:,3)>=spikes(1,3)))/n_rep;
            response(3,id_cluster_SI,cell_n) = spikes(1,3);
            
            %Spatial info analysis
            subplot(4,4,16);    hi = histogram(info,'Normalization','pdf'); hold on;
            stem(info(1),1);   hold off;    title('Spatial Info');
            S_Info(1,id_cluster_SI,cell_n) = info(1);
            S_Info(2,id_cluster_SI,cell_n) = length(find(info>=info(1)))/n_rep;
            
            drawnow();      saveas(gcf,[figures_directory1, '/', batName '_' dateSesh '_' sessionType '_ROI_' num2str(cell_n) '_cluster_' num2str(id_cluster_SI) '.jpg']);
        end
    end
    close all;
end

%% Look at significant cells

if p_val_analysis
    
    %Putative place cells(**during activity & **spatial info & Peak Firing > 3 average firing)
    %pp_cells is matrix of 0 and 1 to show if cell within each cluster
    %(ncluster x ncells)
    pp_cells = [false(1,N); squeeze(S_Info(2,2:end,:)<alfa) & squeeze(p_val(2,2:end,:)<alfa) & squeeze(max(sp_bnd_response(:,2:end,:),[],1))>3*mean(Rate(:,:),2)']; %false for cluster #1
    pp_cells_activity = round(normalize(sp_bnd_response(:,pp_cells),1,'range'),5)'; %normalize activity between 0-1 for putative place fields along 30 space fields
    [sorted_rows,~] = find(pp_cells_activity==1);   sorted_pp_fields = pp_cells_activity(sorted_rows,:);
    figure();   imagesc(sorted_pp_fields,[0 1]);         colormap(viridis); %plots the normalized activity sorted of putative place fields (cells sig/trajectory)
    saveas(gcf,[pwd, '/', batName '_' dateSesh '_' sessionType '_pfields_sorted.jpg']);
    
    %Visualize place fields and calculate centroids
    %outputs every place field as a centroid heatmap along each trajectory
    figure();   set(gcf, 'units','normalized','outerposition',[0.2 0.3 0.5 0.45]);
    [clus,cellNum] = find(pp_cells); %outputs which cells and which clusters have spatial selectivity
    Place_field = struct([]);
    
    %Putative pre cells(**pre activity & Peak Firing > 3 average firing)
    %pp_cells is matrix of 0 and 1 to show if cell within each cluster
    %(ncluster x ncells)
    ppre_cells = [false(1,N); squeeze(p_val(1,2:end,:)<alfa) & squeeze(max(sp_bnd_response_pre(:,2:end,:),[],1))>3*mean(Rate(:,:),2)']; %false for cluster #1
    %ppre_cells_activity = round(normalize(bnd_act_pre(:,ppre_cells),1,'range'),5)'; %normalize activity between 0-1 for putative place fields along 30 space fields
    %[sorted_rows_pre,~] = find(ppre_cells_activity==1);   sorted_ppre_fields = ppre_cells_activity(sorted_rows_pre,:);
    %figure();   imagesc(sorted_ppre_fields,[0 1]);         colormap(viridis); %plots the normalized activity sorted of putative place fields (cells sig/trajectory)
    %saveas(gcf,[pwd, '/', batName '_' dateSesh '_' sessionType '_prefields_sorted.jpg']);
    [clus_pre,cellNum_pre] = find(ppre_cells); %outputs which cells and which clusters have spatial selectivity
    
    %Putative post cells(**pre activity & Peak Firing > 3 average firing)
    %pp_cells is matrix of 0 and 1 to show if cell within each cluster
    %(ncluster x ncells)
    ppost_cells = [false(1,N); squeeze(p_val(3,2:end,:)<alfa) & squeeze(max(sp_bnd_response_pst(:,2:end,:),[],1))>3*mean(Rate(:,:),2)']; %false for cluster #1
    %ppost_cells_activity = round(normalize(bnd_act_pst(:,ppost_cells),1,'range'),5)'; %normalize activity between 0-1 for putative place fields along 30 space fields
    %[sorted_rows_post,~] = find(ppost_cells_activity==1);   sorted_ppost_fields = ppost_cells_activity(sorted_rows_post,:);
    %figure();   imagesc(sorted_ppost_fields,[0 1]);         colormap(viridis); %plots the normalized activity sorted of putative place fields (cells sig/trajectory)
    %saveas(gcf,[pwd, '/', batName '_' dateSesh '_' sessionType '_postfields_sorted.jpg']);
    [clus_post,cellNum_post] = find(ppost_cells); %outputs which cells and which clusters have spatial selectivity
     
    %putative place cells but not looking at spatial info, only based on p-val
    pp_cells_loose = [false(1,N);  squeeze(p_val(2,2:end,:)<alfa) & squeeze(max(sp_bnd_response(:,2:end,:),[],1))>3*mean(Rate(:,:),2)']; %false for cluster #1
    [clus_loose,cellNum_loose] = find(pp_cells_loose);
    
    %putative place cells but not looking at spatial info, only based on
    %p-val and no minimum spike rate
    pp_cells_looseloose = [false(1,N);  squeeze(p_val(2,2:end,:)<alfa)]; %false for cluster #1
    [clus_looseloose,cellNum_looseloose] = find(pp_cells_looseloose);
    
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
        sgtitle(['ROI: ' num2str(cellNum(clus_i)) ' Cluster: ' num2str(clus(clus_i))]);
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
        
        saveas(gcf,[figures_directory1, '/',batName '_' dateSesh '_' sessionType '_Place_field' num2str(clus_i) '.jpg']);
    end
    close all;
    
    %Calculate percentage of place cells
    perc_place = length(unique(cellNum))./N;
    perc_pre = length(unique(cellNum_pre))./N;
    perc_post = length(unique(cellNum_post))./N;
    perc_place_loose = length(unique(cellNum_loose))./N;
    perc_place_looseloose = length(unique(cellNum_looseloose))./N;
    cellNum_none = [unique(cellNum)' unique(cellNum_pre)' unique(cellNum_post)'];
    perc_none = (N-length(unique(cellNum_none)))./N;
    %Calculate percentage of pre, during, post cells on any trajectory (excluding cluster 1)
    perc_any = sum(squeeze(any(p_val(:,2:end,:)<alfa,2)),2)./N;
    
    
    %Evaluate distances between centroids for corresponding place fields, within a
    %single cell
%     field_dist = [];
%     for cell_i = 1:N
%         if length(find([Place_field.cell] == cell_i))>1
%             cell_pairs = nchoosek(find([Place_field.cell] == cell_i),2);
%             for n = 1:size(cell_pairs,1)
%                 field_dist = [field_dist; norm([Place_field(cell_pairs(n,1)).pos]-[Place_field(cell_pairs(n,2)).pos])];
%             end
%         end
%     end
    

    %Look at the activity of significantly modulated cells
    %this may have to be corrected because it is not looking at all
    %pre/post/during cells, just looking at during right now
    pt_cells = [false(1,N); squeeze(any(p_val(:,2:end,:)<alfa))]; %pt cells have signif pval in any period (pre,during, post)
    pt_cells_activity = round(normalize(avg_bnd_act(:,pp_cells),1,'range'),5)';
    [sorted_rows,~] = find(pt_cells_activity==1);   sorted_cells = pt_cells_activity(sorted_rows,:);
    figure();       set(gcf, 'units','normalized','outerposition',[0.25 0.25 0.5 0.5]);
    ax1 = subplot(1,3,1);       imagesc(sorted_cells(:,1:n_bins),[0 1]);         colormap(viridis);
    ax2 = subplot(1,3,2);       imagesc(sorted_cells(:,n_bins+1:2*n_bins),[0 1]);
    ax3 = subplot(1,3,3);       imagesc(sorted_cells(:,2*n_bins+1:end),[0 1]);
    ax1.Title.String = 'Pre-Flight';        ax2.Title.String = 'During-Flight'; ax3.Title.String = 'Post-Flight';
    ax1.XLabel.String = 'Bin';              ax2.XLabel.String = 'Bin';          ax3.XLabel.String = 'Bin';
    ax1.YLabel.String = 'Neuron x Flight #';ax2.YTickLabel= [];                 ax3.YTickLabel= [];
    saveas(gcf,[pwd, '/', batName '_' dateSesh '_' sessionType '_activity_sorted.jpg']);
    
    
    %save the variables in a structure
    placeCells.p_val = p_val;
    placeCells.perc_place = perc_place;
    placeCells.perc_pre = perc_pre;
    placeCells.perc_post = perc_post;
    placeCells.perc_none = perc_none;
    placeCells.perc_place_loose = perc_place_loose;
    placeCells.perc_place_looseloose = perc_place_looseloose;
    placeCells.perc_any = perc_any;
    placeCells.pt_cells = pt_cells;
    placeCells.pt_cells_activity = pt_cells_activity;
    placeCells.pp_cells = pp_cells;
    placeCells.pp_cells_activity = pp_cells_activity;
    placeCells.ppre_cells = ppre_cells;
    placeCells.ppost_cells = ppost_cells;
    placeCells.pp_cells_loose = pp_cells_loose;
    placeCells.pp_cells_looseloose = pp_cells_looseloose;
    placeCells.bnd_act = bnd_act;
    placeCells.bnd_act_pre = bnd_act_pre;
    placeCells.bnd_act_dur = bnd_act_dur;
    placeCells.bnd_act_post = bnd_act_pst;
    placeCells.S_Info = S_Info;
    placeCells.response = response;
    placeCells.avg_bnd_act = avg_bnd_act;
    placeCells.sp_bnd_response = sp_bnd_response;
    placeCells.sp_bnd_velCel = sp_bnd_velCel;
    placeCells.N_cells_total = N;
    placeCells.cellNum_place = cellNum;
    placeCells.cellNum_pre = cellNum_pre;
    placeCells.cellNum_post = cellNum_post;
    placeCells.cellNum_none = cellNum_none;
    placeCells.cellNum_place_loose = cellNum_loose;
    placeCells.cellNum_place_looseloose = cellNum_looseloose;
    placeCells.n_bins = n_bins;
    placeCells.n_rep = n_rep;
    placeCells.alfa = alfa;
    placeCells.pre_dur = pre_dur;
    placeCells.pst_dur = pst_dur;
    placeCells.n_space_bins = n_space_bins;
    placeCells.Ncells = N;
    placeCells.SizeNeuralTrace = T;
    
end

%% Save (modify this to save only relevant variables)
if save_data == 1
    a_filename = [pwd, filesep batName '_' dateSesh '_' sessionType '_ExtractedPlaceCells_',datestr(now, 'yymmdd_HHMM'),'.mat'];
    save(a_filename,'placeCells');
%     if analyze_Ca
%         save([pwd,'/A_', batdate, '.mat'],'A');      %save spatial footprints
%     end
end
disp(datestr(now));