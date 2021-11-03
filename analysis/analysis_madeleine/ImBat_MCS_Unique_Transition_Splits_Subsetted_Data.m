function [TG_transition_number_list,DG_transition_number_list] = ImBat_MCS_Unique_Transition_Splits_Subsetted_Data(CombinedROI,ROI_Data,flightPaths,c_s_34,FlightAlignedROI,cell_set,day_set,cluster_number)

% For a given cell, 
cells2use = cell_set;
%[1, 2, 6, 13, 17, 20, 25, 26, 27, 37, 38, 41, 43, 47, 48, 67, 69, 93]; % ROIs to use Ge 19 to 24
days2use = day_set; % days to use
clu = cluster_number; % cluster (INDEXED BY clu+1)!! (To get cluster 2, clu = 1)

% Plot heatmap alignment of that given cell for those days for that cluster
%ImBat_PlotAlignedROIs(FlightAlignedROI{clu},'cells2use',cells2use);

%% Split those flights into pairs, triplets, or quadruplets
% FG = [];
% for i=4:length(c_s_34)
%     if c_s_34(i) == clu+1 & c_s_34(i-1) == clu+1 & c_s_34(i-2) == clu+1
%         FG = [FG; c_s_34(i-3),c_s_34(i-2),c_s_34(i-1),c_s_34(i)];
%     elseif c_s_34(i-1) == clu+1 
%         FG = [FG; NaN,NaN,NaN,NaN];
%     end
% end
% 
% unique_FG = unique(FG(~isnan(FG)));
% FG_split = cell(length(unique_FG),1);
% FG_split_idx = cell(length(unique_FG),1);
% for i=1:length(unique_FG)
%     for j=1:length(FG)
%         if FG(j,1) == unique_FG(i)
%             FG_split{i} = [FG_split{i};FG(j,:)];
%             FG_split_idx{i} = [FG_split_idx{i};j];
%         end
%     end
% end

% Get all cluster triplets
TG = {};
for j=1:length(unique(c_s_34))
    TG_j = [];
    for i=3:length(c_s_34)
        if c_s_34(i) == clu+1 & c_s_34(i-1) == j
            TG_j = [TG_j; c_s_34(i-2),c_s_34(i-1),c_s_34(i)];
        elseif c_s_34(i-1) == clu+1 
            TG_j = [TG_j; NaN,NaN,NaN];
        end
    end
    TG{j} = TG_j;
end

% Split each TG cell into respective transitions from t-1
for k=1:length(unique(c_s_34))
    t_tg = TG{k};
    unique_TG = unique(t_tg(~isnan(t_tg)));
    TG_split = cell(length(unique_TG),1);
    TG_split_idx = cell(length(unique_TG),1);
    for i=1:length(unique_TG)
        for j=1:length(t_tg)
            if t_tg(j,1) == unique_TG(i)
                TG_split{i} = [TG_split{i};t_tg(j,:)];
                TG_split_idx{i} = [TG_split_idx{i};j];
            end
        end
    end
    TG_split_agg{k} = TG_split;
    TG_split_idx_agg{k} = TG_split_idx; 
end
% Each transition type 
TG_transition_number_list = [];
for i=1:length(TG_split_idx_agg)
    tt_temp = cellfun(@numel,TG_split_agg{i})./3;
    TG_transition_number_list = [TG_transition_number_list;tt_temp];
end
% Check
sum(TG_transition_number_list) == length(c_s_34(c_s_34 == clu+1))

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
DG_transition_number_list = cellfun(@numel,DG_split)./2;

end