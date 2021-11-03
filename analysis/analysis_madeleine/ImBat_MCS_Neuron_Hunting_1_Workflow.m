%% Load in required structs (c_s_34 and flightpaths)

% Concatonate/Cluster data across days
[CombinedROI,ROI_Data] = ImBat_GroupCalcium(ROI_Data,cell_registered_struct,aligned_data_struct);

% repair bat flight data ( important for z bats with bad tracking) :
ROI_Data = ImBat_RepairFlightData(ROI_Data);

% now, we have a restricted set for ROI_Data, now cluster the flights:
flightPaths = ImBat_GroupFlights(ROI_Data,'mtf',master_track_file,'dist',1.1);         % just the flights

% Sub in the c_s_34 overall clusters for the subset clusters
[ss_1,rr_1] = sort(flightPaths.flight_starts_idx);
c_s_34 = flightPaths.id(rr_1);

%% Align all flights to the top 10 flight clusters 

for flight_cluster = 1:17
    [FlightAlignedROI{flight_cluster}] = ImBat_Align_FC(CombinedROI,flightPaths,flight_cluster+1);
end 

%% Look at a given neuron's response over all flights of a given cluster for a set of days

cells2use = [1, 2, 6, 13, 17, 20, 25, 26, 27, 37, 38, 41, 43, 47, 48, 67, 69, 93]; % ROIs to use Ge 19 to 24
days2use = 1:5; % days to use
clu = 2; % cluster (INDEXED BY clu+1)!! (To get cluster 2, clu = 1)
ImBat_PlotAlignedROIs(FlightAlignedROI{clu},'cells2use',cells2use);

%% Split those flights into pairs, triplets, or quadruplets
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

%% Run ImBat_MCS_Neuron_Counting to get the timeline of when those preceeding flights happen

%% Plot the splits of DG, TG, or FG
clu_split=5; % This is the desired group of flights (222, or 422, etc)
clear subset Flight_Aligned_ROIs_subset subset_idx
subset = ones(length(DG),1);
subset(DG_split_idx{clu_split}) = 0;
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

%% Look at the PST to see if there is a chunk you should be looking at 
[out_markov] = ImBat_New_Markov(flightPaths);
ImBat_ProbSuffixTree(out_markov,6);

% Stats for above ^
[Out_Markov] = ImBat_PlotMarkov(out_markov,FlightAlignedROI{clu},cells2use);
[p_val p_val_pre] = ImBat_HistoryEncode(Out_Markov);




            

