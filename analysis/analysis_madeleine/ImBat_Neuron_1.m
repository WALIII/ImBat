% Look at a stretch of days look at the cells

% For this stretch of days, find the overall clustering c_s_34 vector 
%load('../Processed_Ge_LongHaul/c_s_34.mat');
%load('../Processed_Ge_LongHaul/flightPaths34.mat');
%[ss,rr] = sort(flightPaths34.flight_starts_idx);
%day_idx = flightPaths34.day(rr);

%load('aligned_data_struct.mat');
%load('cellRegistered_20210126_202011.mat');
%load('ROI_Data.mat');

date_list = cell2mat(flightPaths34.Dates');
for i=1:length(date_list)
    date_nums(i) = str2double(date_list(i,:));
end

for i=1:length(ROI_Data)
    date_index_list(i) = find(date_nums == str2double(ROI_Data{i}.date(3:end)))
end

% Stretch of c_s_34
c_s_34_seg_idxs = [];
for i=1:length(date_index_list)
    temp = find(day_idx == date_index_list(i))
    c_s_34_seg_idxs = [c_s_34_seg_idxs; c_s_34(temp)];
end
  

% Load in the neurons!
% Concatonate/Cluster data across days
[CombinedROI,ROI_Data] = ImBat_GroupCalcium(ROI_Data,cell_registered_struct,aligned_data_struct);

% repair bat flight data ( important for z bats with bad tracking) :
ROI_Data = ImBat_RepairFlightData(ROI_Data);

% now, we have a restricted set for ROI_Data, now cluster the flights:
flightPaths = ImBat_GroupFlights(ROI_Data,'mtf',master_track_file,'dist',1.1);         % just the flights

% Sub in the c_s_34 overall clusters for the subset clusters
[ss_1,rr_1] = sort(flightPaths.flight_starts_idx);
c_s_34 = flightPaths.id(rr_1);
%flightPaths.id = ordered_flight_ids;
%group the cluster ids
%for i = 1:max(flightPaths.id)
%    flightPaths.clusterIndex{i} = find(flightPaths.id == i);
%end
    
% labeled dff projections:
days2use = 1:5; % days to use
cells2use = 41;%[1, 2, 6, 13, 17, 20, 25, 26, 27, 37, 38, 41, 43, 47, 48, 67, 69, 93]; % ROIs to use Ge 19 to 24
ImBat_PlotTrackedMasks(ROI_Data,CombinedROI,cell_registered_struct,days2use,cells2use);

% Align Flight data to the top 3 flight clusters:
for flight_cluster = 1:4
    [FlightAlignedROI{flight_cluster}] = ImBat_Align_FC(CombinedROI,flightPaths,flight_cluster+1);
end 
  
%% For flight cluster clu, plot the ROIs aligned to takeoff
clu = 1;
ImBat_PlotAlignedROIs(FlightAlignedROI{clu},'cells2use',cells2use);

%% Get all the cluster pairs for clu
FG = [];
for i=4:length(c_s_34)
    if c_s_34(i) == clu+1 & c_s_34(i-1) == clu+1 & c_s_34(i-2) == clu+1
        FG = [FG; c_s_34(i-3),c_s_34(i-2),c_s_34(i-1),c_s_34(i)];
    elseif c_s_34(i-1) == clu+1 
        FG = [FG; NaN,NaN,NaN,NaN];
    end
end

unique_FG = unique(FG(~isnan(FG)));
FG_split = cell(length(unique_FG),1);
FG_split_idx = cell(length(unique_FG),1);
for i=1:length(unique_FG)
    for j=1:length(FG)
        if FG(j,1) == unique_FG(i)
            FG_split{i} = [FG_split{i};FG(j,:)];
            FG_split_idx{i} = [FG_split_idx{i};j];
        end
    end
end

TG = []; %TG = TGM_pruned_3;
for i=3:length(c_s_34)
    if c_s_34(i) == clu+1 & c_s_34(i-1) == clu+1
        TG = [TG; c_s_34(i-2),c_s_34(i-1),c_s_34(i)];
    elseif c_s_34(i-1) == clu+1 
        TG = [TG; NaN,NaN,NaN];
    end
end

unique_TG = unique(TG(~isnan(TG)));
TG_split = cell(length(unique_TG),1);
TG_split_idx = cell(length(unique_TG),1);
for i=1:length(unique_TG)
    for j=1:length(TG)
        if TG(j,1) == unique_TG(i)
            TG_split{i} = [TG_split{i};TG(j,:)];
            TG_split_idx{i} = [TG_split_idx{i};j];
        end
    end
end

% Get all the cluster pairs
DG = [];
for i=3:length(c_s_34)
    if c_s_34(i) == clu+1
        DG = [DG; c_s_34(i-1),c_s_34(i)];
    end
end

unique_DG = unique(DG);
DG_split = cell(length(unique_DG),1);
DG_split_idx = cell(length(unique_DG),1);
for i=1:length(unique(DG))
    for j=1:length(DG)
        if DG(j,1) == unique_DG(i)
            DG_split{i} = [DG_split{i};DG(j,:)];
            DG_split_idx{i} = [DG_split_idx{i};j];
        end
    end
end

% Make a kernal smoothing plot of when these things occur to look for when to take out
% Kernal smooth plot
figure(); hold on;
colors_1 = jet(length(unique_TG));
for i=2:4%length(unique(TG))
    pdSix = fitdist(TG_split_idx{i},'Kernel','BandWidth',4);
    x = 0:.1:45;
    ySix = pdf(pdSix,x);
    plot(x,ySix,'k-','LineWidth',2,'Color',colors_1(i,:));
    legend();
end

% Create seperate FlightROI Structure and plot whatever preceding cluster
% subset you want. In this case, I'm subsetting the cluster 2 flights into 
% 222, 422, 522
clu_split=4;
clear subset Flight_Aligned_ROIs_subset subset_idx
subset = ones(length(FG),1);
subset(FG_split_idx{clu_split}) = 0;
subset_idx = find(subset==1);
Flight_Aligned_ROIs_subset = FlightAlignedROI{clu};
Flight_Aligned_ROIs_subset.C(:,:,subset_idx) = [];
Flight_Aligned_ROIs_subset.C_raw(:,:,subset_idx) = [];
Flight_Aligned_ROIs_subset.S(:,:,subset_idx) = [];
Flight_Aligned_ROIs_subset.IDX(subset_idx) = [];
Flight_Aligned_ROIs_subset.cluster_idX(subset_idx) = [];
Flight_Aligned_ROIs_subset.ClustFlight_withPads(:,:,subset_idx) = [];
Flight_Aligned_ROIs_subset.ClustFlight(:,:,subset_idx) = [];
Flight_Aligned_ROIs_subset.FlightLength(subset_idx) = [];
Flight_Aligned_ROIs_subset.FlightTimes(subset_idx) = [];
Flight_Aligned_ROIs_subset.CutCells_date(subset_idx) = [];

ImBat_PlotAlignedROIs(Flight_Aligned_ROIs_subset,'cells2use',cells2use);

            
            
            