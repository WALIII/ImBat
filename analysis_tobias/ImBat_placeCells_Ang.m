function [plotFiringTrajectory] = ImBat_PlaceCells_Tobias(flightPaths, cellData, alignment,varargin)

batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
p_val_analysis = 1;

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
n_bins = 10;     %good around here                                                           %number of bins to divide the pre-during-post flight interval
n_rep = 1000;    %can lower to 10 for debugging                                                           %number of shufflings
alfa = 0.05;     %can lower to 0.01 or 0.1                                                           %significance level
pre_dur = 3;     %play with 3-5                                                           %duration of the pre flight period (s):     comparable with flight dur
pst_dur = 3;     %play with 3-5                                                           %duration of the post flight period (s):    but shorter than half interflight
w = gausswin(1); %keep at 1 for no smoothing                                                           %witdh of the gaussian filter (number of bins), use odd values. 1=no filtering
n_space_bins = 30;  %correlates with space resolution but good around 30 (20cm chunks)                                                        %number of spatial bins

%Flight Room references (provvisory)
xR = +2.85; xL = -2.85; yF = 2.50;  yB = -2.50;  zT = 2.20;                 %Flight volume coordinates
F1 = [2.56; 1.23; 1.72];    F2 = [2.56; -1.04; 1.72];                       %Feeders coordinates
F3 = [2.56; 1.43; 0.72];    F4 = [2.56; -1.24; 0.72];                       %Feeders coordinates
edges_d = {xL:(xR-xL)/10:xR yB:(yF-yB)/10:yF};                              %Edges for density histogram

%Create analysis folder for storing the results
analysis_directory=fullfile(pwd,['Analysis_',datestr(now, 'yymmdd_HHMM')]);
if ~exist(analysis_directory,'dir')
    mkdir(analysis_directory);
end
%% Code in the middle



%% --------------------------------- REFINED ANALYSIS START FROM HERE------------------------------
%//////////////////////////////////////////////////////////////////////////////////////////////////
% Calculations are done on the Rate matrix, which is the S matrix,
% normalized by standard deviation and smoothed on a 0.5s window

if p_val_analysis
    until_cluster = min(n_surv_clusters,100);   %change this if you want less clusters
    
    %Create folder to store figures
    figures_directory1=fullfile(analysis_directory,'Spatial_information');
    delete([figures_directory1 '\*']);
    if exist(figures_directory1,'dir')~=7
        mkdir(figures_directory1);
    end
    
    %Initialization of matrices and arrays
    frames_to_shift(1)=0;     frames_to_shift(2:n_rep) = randi([10*CNMF_Fs T-10*CNMF_Fs],1,n_rep-1);   %Shifting in time (longer than 10s)
    p_val = zeros(3,until_cluster,N); %pre(1),during(2),post(3)/#cluster/#nueron                       %p values for pre, during, post active neurons
    response = zeros(3,until_cluster,N); %sum of each phase                   %Integrated response during pre, during, post periods 'sum(median(Rate)'
    avg_bnd_act = zeros(3*n_bins,until_cluster,N);   %n_bins per section, for each pre/dur/post and for each cluster & cell       %Activity across bins from pre to post
    sp_bnd_response = zeros(n_space_bins,until_cluster,N);  %Spatially binned activity along the trajectory
    S_Info = zeros(2,until_cluster,N); %2 b/c 1dim=actual info, 2dim=p-val                     %bits and p value for spatial information
    
    %Binning in time and space, p values calculation
    figure();   set(gcf, 'units','normalized','outerposition',[0.2 0 0.5 1]);
    for id_cluster_SI = 1:until_cluster %for each cluster
        id = [];    id = find(flight_clus.id==id_cluster_SI); %find all flights that belong to that cluster
        for cell_n = 1:N %for each cell, initialize the below matrices
            sgtitle(['ROI: ' num2str(cell_n) ' Cluster: ' num2str(id_cluster_SI)]);
            disp(['Cell number: ' num2str(cell_n) '  Trajectory number: ' num2str(id_cluster_SI)]);
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
                
                for ii=1:size(id,1) %for all flights within the cluster, define the following vectors
                    %grabbing and binning the velocity and activity along
                    %with other variables
                    Act_pre = [];   Act_dur = [];   Act_pst = [];   v_trj = [];
                    Act_pre =  Rate_sh(flight_clus_ds.strt_frame(id(ii))-pre_dur*CNMF_Fs:flight_clus_ds.strt_frame(id(ii))-1);
                    Act_dur =  Rate_sh(flight_clus_ds.strt_frame(id(ii)):flight_clus_ds.stop_frame(id(ii)));
                    Act_pst =  Rate_sh(flight_clus_ds.stop_frame(id(ii))+1:flight_clus_ds.stop_frame(id(ii))+pst_dur*CNMF_Fs);
                    
                    %Temporally binned activity for pre-during-post
                    flight_dur(1,ii) = flight_clus_ds.dur(id(ii)); %this comes from the flightPaths output
                    bnd_act_pre(:,ii) = interp1(linspace(1,100,size(Act_pre,1)),Act_pre,linspace(1,100,n_bins),'linear')';
                    bnd_act_dur(:,ii) = interp1(linspace(1,100,size(Act_dur,1)),Act_dur,linspace(1,100,n_bins),'linear')';
                    bnd_act_pst(:,ii) = interp1(linspace(1,100,size(Act_pst,1)),Act_pst,linspace(1,100,n_bins),'linear')';
                    
                    %Spatially binned activity for spatial information
                    %during flight (spatial activity) and interpolating
                    %across the bins
                    v_trj = flight_clus_ds.vel(1,~isnan(flight_clus_ds.pos(1,:,id(ii))),id(ii));
                    sp_bnd_act(:,ii) = interp1(linspace(1,flight_clus_ds.length(id(ii)),size(Act_dur,1)),Act_dur,linspace(1,flight_clus_ds.length(id(ii)),n_space_bins),'linear')';
                    sp_bnd_vel(:,ii) = interp1(linspace(1,flight_clus_ds.length(id(ii)),size(v_trj,2)),v_trj',linspace(1,flight_clus_ds.length(id(ii)),n_space_bins),'linear')';
                    sp_bnd_path(:,ii) = linspace(1,flight_clus_ds.length(id(ii)),n_space_bins)';
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
                    avg_bnd_act(:,id_cluster_SI, cell_n) = filter(w,1,median(bnd_act,2));
                    sp_bnd_response(:,id_cluster_SI,cell_n) = filter(w,1,lambda);
                    
                    ciplot(filter(w,1,median(bnd_act,2))-std(bnd_act,[],2)./sqrt(size(id,1)),filter(w,1,median(bnd_act,2))+std(bnd_act,[],2)./sqrt(size(id,1)),linspace(1,100,3*n_bins));
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
                spikes(n,1) = sum(median(bnd_act_pre,2))/pre_dur;
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
            
            drawnow();      saveas(gcf,[figures_directory1, '/', 'ROI_' num2str(cell_n) '_cluster_' num2str(id_cluster_SI) '.png']);
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
    saveas(gcf,[analysis_directory, '/', batdate '_pfields_sorted.png']);
    
    %Visualize place fields and calculate centroids
    %outputs every place field as a centroid heatmap along each trajectory
    figure();   set(gcf, 'units','normalized','outerposition',[0.2 0.3 0.5 0.45]);
    [clus,cell] = find(pp_cells); %outputs which cells and which clusters have spatial selectivity
    Place_field = struct([]);
    for i = 1:length(clus) %for each cluster
        
        %take each flight from a place cell cluster, take the shortest
        %flight, and generate an average trajectory based off the shortest
        %flight
        id = find(flight_clus_ds.id==clus(i));
        shortest_flight_i = min(squeeze(sum(~isnan(flight_clus_ds.pos(:,:,id)),2)),[],2);
        ave_trajectory = nanmean(flight_clus_ds.pos(:,1:shortest_flight_i(1),id),3);
        %ave_acceleratn = nanmean(flight_clus_ds.acc(1,1:shortest_flight_i(1),id),3);
        
        npo = size(ave_trajectory,2); %length of average trajectory
        take_off = mean(squeeze(flight_clus.pos(:,1,id)),2);
        map = interp1(linspace(1,100,size(sp_bnd_response(:,clus(i),cell(i)),1)),sp_bnd_response(:,clus(i),cell(i)),linspace(1,100,npo),'linear')';
        map_color = uint8(round(normalize(map,'range',[2 100]))); %make the mapping between activity and color
        
        %Determine location of the max activity 
        [~,max_bin] = max(sp_bnd_response(:,clus(i),cell(i)),[],1,'linear');
        max_position = round(size(ave_trajectory,2)*max_bin/n_space_bins); %find centroid of place field
        %plot the firing activity along the average 
        sgtitle(['ROI: ' num2str(cell(i)) ' Cluster: ' num2str(clus(i))]);
        subplot(121); cmap = viridis(100);
        p = plot3(ave_trajectory(1,:),ave_trajectory(2,:),ave_trajectory(3,:),'r', 'LineWidth',5);
        grid on;
        cdata = [uint8(cmap(map_color,:)*255) uint8(ones(npo,1))].';
        drawnow();  set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cdata);
        hold on;
        textscatter3(take_off(1),take_off(2),take_off(3),"Take-off");
        scatter3(ave_trajectory(1,max_position),ave_trajectory(2,max_position),ave_trajectory(3,max_position));
        xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
        xlabel('x');    ylabel('y');    zlabel('z');
        drawnow();  hold off;
        
        subplot(122);
        plot([1:1:n_space_bins],sp_bnd_response(:,clus(i),cell(i)),'LineWidth',3);   xlabel('Space along trajectory');    ylabel('Activity (SD units)');
        xticks([1 n_space_bins]);   xticklabels({'Take-off','Landing'});
        %save info for each place field (position, cell num, clust num)
        Place_field(i).pos = ave_trajectory(:,max_position);
        Place_field(i).cell = cell(i);
        Place_field(i).clus = clus(i);
        
        saveas(gcf,[figures_directory, '/', 'Place_field' num2str(i) '.png']);
    end
    close all;
    
    %Calculate percentage of place cells
    perc_place = length(unique(cell))./N;
    
    %Evaluate distances between centroids for corresponding place fields, within a
    %single cell
    field_dist = [];
    for i = 1:N
        if length(find([Place_field.cell] == i))>1
            cell_pairs = nchoosek(find([Place_field.cell] == i),2);
            for n = 1:size(cell_pairs,1)
                field_dist = [field_dist; norm([Place_field(cell_pairs(n,1)).pos]-[Place_field(cell_pairs(n,2)).pos])];
            end
        end
    end
    
    %Calculate percentage of pre, during, post cells on any trajectory (excluding cluster 1)
    percentage = sum(squeeze(any(p_val(:,2:end,:)<alfa,2)),2)./N;
    
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
    saveas(gcf,[analysis_directory, '/', batdate '_activity_sorted.png']);
    
end

%% Save (modify this to save only relevant variables)
if save_data
    a_filename = [analysis_directory,'/Extracted_trajectories_&_activity',datestr(now, 'yymmdd_HHMM'),'.mat'];
    save(a_filename);
    if analyze_Ca
        save([analysis_directory,'/A_', batdate, '.mat'],'A');      %save spatial footprints
    end
end
disp(datestr(now));