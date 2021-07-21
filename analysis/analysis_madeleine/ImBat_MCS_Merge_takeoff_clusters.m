function [KC_new,unique_takeoffs_new] = ImBat_MCS_Merge_takeoff_clusters(KC,unique_takeoffs,user_input);

% Merge the clusters specified in "clusters_to_merge" and
% "dominant_cluster"
if contains(user_input,"merge") & ~contains(user_input,"split")
    
    dominant_cluster = str2double(user_input(end));
    uip = user_input(7:end-2);
    for i=1:size(uip,2)
        clusters_to_merge(i) = str2double(uip(i));
    end
    clusters_to_merge = clusters_to_merge(~isnan(clusters_to_merge));

    % Create a new KC by replacing all of the clusters_to_merge with
    % dominant_cluster
    KC_new = KC;
    KC_replaced_idx = zeros(size(KC,1),1);
    for i=1:size(KC,1)
        if ismember(KC(i),clusters_to_merge)
            KC_replaced_idx(i) = 1;
            KC_new(i) = dominant_cluster;
        end
    end

    % Create a new unique_takeoffs
    unique_takeoffs_new = unique_takeoffs;
    for i=2:size(unique_takeoffs_new,2)
        if size(cell2mat(unique_takeoffs_new(i)),1) > 1
            unique_takeoffs_new(i) = {dominant_cluster};
        end
    end
    
elseif contains(user_input,"merge") & contains(user_input,"split")
    disp("Build");
elseif ~contains(user_input,"merge") & contains(user_input,"split")
    disp("Build");
end



end