function [co,number_unique_transitions_2,number_unique_transitions_3] = ImBat_MCS_Count_Unique_Transitions(c_s_34,flightPaths34,minimum_number_of_transition_instances)

% Determine the flight index cutoff that comprises 95% of the data
c_s_34_seg = c_s_34;
five_percent = round(size(c_s_34,1)*(0.05));
for i=1:length(unique(c_s_34))-1
    temp1 = length(c_s_34(c_s_34 > i));
    temp2 = length(c_s_34(c_s_34 > (i+1)));
    if temp2 < five_percent & temp1 > five_percent
        co = i;
        break
    else
        co = 0;
    end
end     

%% Extract every 2-gram possible from the data using sliding window
DGM_all = [];
window_size=1;
for i=1:length(c_s_34_seg)-window_size
    two_gram = c_s_34_seg(i:i+window_size);
    DGM_all(i,:) = two_gram;
end

% Remove all 2-grams that bridge sessions
DGM_blacklist = [];
[ss,rr] = sort(flightPaths34.flight_starts_idx);
day_idx = flightPaths34.day(rr);
dur_idx = flightPaths34.dur(rr);
for i=1:length(DGM_all)-2
    if day_idx(i) ~= day_idx(i+1)
        DGM_blacklist = [DGM_blacklist,[i-1,i,i+1]];
    end
end
DGM_pruned_1 = DGM_all;
DGM_pruned_1(DGM_blacklist,:) = [];

% Remove all 2-grams that contain flights greater than the cutoff 
DGM_pruned_2 = [];
for i=1:length(DGM_pruned_1)
    if sum(DGM_pruned_1(i,:) > co) >= 1
        continue
    else
        DGM_pruned_2 = [DGM_pruned_2;DGM_pruned_1(i,:)];
    end
end

% Remove all 3-grams containing a 1
DGM_pruned_3 = [];
for i=1:length(DGM_pruned_2)
    if sum(DGM_pruned_2(i,:)==1) >= 1 
        continue
    else
        DGM_pruned_3 = [DGM_pruned_3;DGM_pruned_2(i,:)];
    end
end

DGM_T_Mat = zeros(length(unique(DGM_pruned_3)),length(unique(DGM_pruned_3)));
for i=1:length(unique(DGM_pruned_3))
    for j=1:length(unique(DGM_pruned_3))
        for m=1:length(DGM_pruned_3)
            if DGM_pruned_3(m,1) == i & DGM_pruned_3(m,2) == j
                DGM_T_Mat(i,j) = DGM_T_Mat(i,j) + 1;
            end
        end
    end
end
%figure(); heatmap(DGM_T_Mat)
DGM_T_Mat_Vec = reshape(DGM_T_Mat,length(DGM_T_Mat)*length(DGM_T_Mat),1);

number_unique_transitions_2 = length(DGM_T_Mat_Vec(DGM_T_Mat_Vec>minimum_number_of_transition_instances));
unique_transitions_values_2 = find(DGM_T_Mat_Vec>minimum_number_of_transition_instances);
%figure();hold on; bar(DGM_T_Mat_Vec); yline(minimum_number_of_transition_instances);
    

%% Extract every 3-gram possible from the data using sliding window
c_s_34_seg = c_s_34
TGM_all = [];
window_size=2;
for i=1:length(c_s_34_seg)-window_size
    three_gram = c_s_34_seg(i:i+window_size);
    TGM_all(i,:) = three_gram;
end

% Calculate co based on which flight types make up >5% of the data
five_percent = round(size(c_s_34,1)*(0.05));
colormap = distinguishable_colors(co);
length(c_s_34(c_s_34 > co))

% Remove all 3-grams that bridge sessions
TGM_blacklist = [];
[ss,rr] = sort(flightPaths34.flight_starts_idx);
day_idx = flightPaths34.day(rr);
dur_idx = flightPaths34.dur(rr);
for i=1:length(TGM_all)-2
    if day_idx(i) ~= day_idx(i+1)
        TGM_blacklist = [TGM_blacklist,[i-1,i,i+1]];
    end
end
TGM_pruned_1 = TGM_all;
TGM_pruned_1(TGM_blacklist,:) = [];

% Remove all 3-grams that contain flights greater than the cutoff 
TGM_pruned_2 = [];
for i=1:length(TGM_pruned_1)
    if sum(TGM_pruned_1(i,:) > co) >= 1
        continue
    else
        TGM_pruned_2 = [TGM_pruned_2;TGM_pruned_1(i,:)];
    end
end

% Remove all 3-grams containing a 1
TGM_pruned_3 = [];
for i=1:length(TGM_pruned_2)
    if sum(TGM_pruned_2(i,:)==1) >= 1 
        continue;
    else
        TGM_pruned_3 = [TGM_pruned_3;TGM_pruned_2(i,:)];
    end
end

% Calculate 3-count from the TGM matrix data
TGM_T_MATRIX = zeros(co*co,co);
for i=1:co  
    for j=1:co 
        for m=1:length(TGM_pruned_3)                          
            for k=1:co                                  
                if TGM_pruned_3(m,3) == j & TGM_pruned_3(m,2)== i & TGM_pruned_3(m,1) == k    
                    TGM_T_MATRIX((i-1)*co + j,k) = TGM_T_MATRIX((i-1)*co + j,k) + 1;                 
                end
            end
        end
    end
end
    
%figure(); heatmap(TGM_T_MATRIX);
TGM_T_MATRIX_Vec = reshape(TGM_T_MATRIX,length(TGM_T_MATRIX)*size(TGM_T_MATRIX,2),1);

number_unique_transitions_3 = length(TGM_T_MATRIX_Vec(TGM_T_MATRIX_Vec>minimum_number_of_transition_instances));
unique_transitions_values_3 = find(TGM_T_MATRIX_Vec>minimum_number_of_transition_instances);
%figure(); hold on; bar(TGM_T_MATRIX_Vec); yline(minimum_number_of_transition_instances); 
 
end