function [takeoff_locations_cpu,pruned_dataset,KC,outliers] = ImBat_MCS_determine_takeoff_locations(flightPaths,c_s,k)

    % Produces 
    %   pruned_dataset.mat - contains the indexes of which flights to use (that are not spatial outliers)
    %   takeoff_locations_cpu.mat - contains the list of takeoff locations that corrospond to the flight types in PRUNED c_s
    % When you get a good clustering, save the takeoff_locations_cpu and
    % pruned_dataset.mat files
    
    
    % Sort the start indexes, xyz positions, and velocities according to
    % the AllFlightsMasterTime vector.
    flight_start_times = flightPaths.AllFlightsMasterTime(flightPaths.flight_starts_idx(:));
    [tts,tts_idx]  = sort(flight_start_times(:));

    flight_starts_xyz_unsorted = squeeze(flightPaths.flight_starts_xyz);
    flight_starts_xyz = flight_starts_xyz_unsorted(:,tts_idx);

%     figure();
%     scatter(flight_starts_xyz_unsorted(1,:),flight_starts_xyz_unsorted(2,:),flight_starts_xyz_unsorted(3,:));
%     figure();
%     scatter(flight_starts_xyz(1,:),flight_starts_xyz(2,:),flight_starts_xyz(3,:));

    flight_pos_xyz = flightPaths.pos(:,:,tts_idx);

    fd = flightPaths.day(tts_idx);

    flight_vel_unsorted = squeeze(flightPaths.vel);
    flight_vel = flight_vel_unsorted(:,tts_idx);

%     figure(); 
%     scatter([1:size(flight_vel_unsorted,2)],flight_vel_unsorted(1,:));
%     figure();
%     scatter([1:size(flight_vel,2)],flight_vel(1,:));

    % Identify the spatially-weird outliers (NOT SORTED BY TIME)
    [outlier_flight_indexes]=ImBat_MCS_find_outlier_flights(flightPaths);
    % For now, exclude these points

    % Identify the velocity-weird outliers
    [outlier_velocity_indexes] = ImBat_MCS_find_velocity_outlier_flights(flightPaths,outlier_flight_indexes);
    % For now, include these points
    
    % Remove the indexes of the outlier flights (from the SORTED flight start
    % xys position vector)
    bin_outlier_flight_indexes = find(outlier_flight_indexes==1);
    non_outlier_indexes= [1:size(tts_idx,1)];
    non_outlier_indexes(bin_outlier_flight_indexes) = [];
    flight_starts_xyz_pruned = flight_starts_xyz;
    flight_starts_xyz_pruned(:,bin_outlier_flight_indexes) = [];
    flight_pos_xyz_pruned = flight_pos_xyz;
    flight_pos_xyz_pruned(:,:,bin_outlier_flight_indexes) = [];

    % Cluster the non-outlier points according to Kmeans
    clear Xvar Yvar Zvar;
    max_clusters=100;
    for j=1:max_clusters
        C = kmeans(flight_starts_xyz_pruned',j);
        for i=1:max(C)
            Xvar(j,i) = var(flight_starts_xyz_pruned(1,find(C==i)));
            Yvar(j,i) = var(flight_starts_xyz_pruned(2,find(C==i)));
            Zvar(j,i) = var(flight_starts_xyz_pruned(3,find(C==i)));
        end
    end
    Xvar(Xvar==0) = NaN;
    Yvar(Yvar==0) = NaN;
    Zvar(Zvar==0) = NaN;
    figure(); hold on;
    title("Mean variance of X, Y, and Z coordinates within a cluster for 1:100 clusters");
    plot([1:max_clusters],nanmean(Xvar,2));
    plot([1:max_clusters],nanmean(Yvar,2));
    plot([1:max_clusters],nanmean(Zvar,2));

    % It seems like 8 clusters is nice for this dataset? 
    colormap = jet(k);
    KC = kmeans(flight_starts_xyz_pruned',k);
    figure();
    hold on;
    colormap = jet(k);
    xlim([-3 3]); ylim([-2.5 2.5]); zlim([0,2.5]);
    for i=1:max(KC)
        scatter3(flight_starts_xyz_pruned(1,find(KC==i)),flight_starts_xyz_pruned(2,find(KC==i)),flight_starts_xyz_pruned(3,find(KC==i)),[], colormap(i,:),'filled','MarkerEdgeColor','k');
    end

   % Check your work silly goose by plotting the starting points of the
   % outlier flights and the rest of those flights in blue
    figure(); hold on; title("Spatial outlier flights. Starting point in red");
    for i=1:size(bin_outlier_flight_indexes,2)
        FN=bin_outlier_flight_indexes(i);
        plot3(flight_pos_xyz(1,:,FN),flight_pos_xyz(2,:,FN),flight_pos_xyz(3,:,FN),'o','Color','b','MarkerSize',2);
    end    
    for i=1:size(bin_outlier_flight_indexes,2)
        FN=bin_outlier_flight_indexes(i);
        plot3(flight_starts_xyz(1,FN),flight_starts_xyz(2,FN),flight_starts_xyz(3,FN),'o','Color','r','MarkerSize',10);
    end 
    
    % Now match the flight type to the takeoff location
    % For every flight type, how many different takeoff locations are there?
    flight_cluster_ids = c_s;%flightPaths34.id;
    flight_cluster_ids(bin_outlier_flight_indexes)=[];
    unique_takeoffs = {};    
    for i=1:size(unique(flight_cluster_ids),1)
        flt_idxs = find(flight_cluster_ids == i);
        tkff_vls = KC(flt_idxs);
        unique_takeoffs{i} = unique(tkff_vls);
    end
    
    % Plot number of possible takeoff location clusters by number of flight
    % types that have multiple potential takeoff locations
    clear takeoffs_per_flight_type multi_takeoff_flights
    max_clusters=100;
    for j=1:max_clusters
        C = kmeans(flight_starts_xyz_pruned',j);
        unique_takeoffs = {};    
        for i=1:size(unique(flight_cluster_ids),1)
            flt_idxs = find(flight_cluster_ids == i);
            tkff_vls = C(flt_idxs);
            unique_takeoffs{i} = unique(tkff_vls);
        end
        temp_vec = cellfun('length',unique_takeoffs)';
        takeoffs_per_flight_type(j,:) = temp_vec;
        multi_takeoff_flights(j) = sum(temp_vec>1)-1;
    end 
    
    figure(); title("Number of flight types with multiple takeoff locations");
    plot(multi_takeoff_flights);
    
    % Determine what flight types the outliers were
    outlier_flight_types = c_s(bin_outlier_flight_indexes);
    disp("Flight Types that are Outliers:");
    disp(outlier_flight_types');
    percent_outliers_areone = sum(outlier_flight_types==1)/size(outlier_flight_types,1);
    disp(strcat(num2str(percent_outliers_areone*100),"% of outliers are cluster 1")); 
    
    % Now figure out which clusters are the starts of which flight types. 
    unique_takeoffs = {};  
    %KC = kmeans(flight_starts_xyz_pruned',k);
    for i=1:size(unique(flight_cluster_ids),1)
        flt_idxs = find(flight_cluster_ids == i);
        tkff_vls = KC(flt_idxs);
        unique_takeoffs{i} = unique(tkff_vls);
    end
    unique_takeoffs{1} = [];
    
    % Plot all the clustered flight takeoffs for clarity
    figure(); hold on; title("Clustered Takeoff Locations, Excluding Spatial Outliers.");
    for i=1:max(KC)
        scatter3(flight_starts_xyz_pruned(1,find(KC==i)),flight_starts_xyz_pruned(2,find(KC==i)),flight_starts_xyz_pruned(3,find(KC==i)),[], colormap(i,:),'filled','MarkerEdgeColor','k');
    end
    
    % See how many points per cluster (color matched)
    for i=1:max(KC)
        pts_per_cluster(i) = sum(KC==i);
    end
    figure(); hold on; b=bar(pts_per_cluster,'FaceColor','flat');
    title("Quantity of flights for each starting location");
    for i=1:max(KC)
        b.CData(i,:) = colormap(i,:);
    end
    
    % Look at the cases where one flight type has two takeoff locations
    for i=1:size(unique_takeoffs,2)
        if size(unique_takeoffs{i},1) > 1
            figure(); hold on;
            % Plot all the starting points of this flight type
            % color-divided by cluster identity
            flt_idxs = find(flight_cluster_ids == i);
            clust_idxs = KC(flt_idxs);
            clusts = unique_takeoffs{i};
            title(strcat("Flight Type ",num2str(i)," has multiple takeoff locations: ",num2str(clusts)));
            cmap = jet(size(clusts,1));
            for j=1:size(unique_takeoffs{i},1)
                FN = flt_idxs(find(clust_idxs == clusts(j)));
                for k=1:size(FN,1)
                    plot3(flight_pos_xyz_pruned(1,:,FN(k)),flight_pos_xyz_pruned(2,:,FN(k)),flight_pos_xyz_pruned(3,:,FN(k)),'o','Color',[0.7 0.7 0.7],'MarkerSize',3);
                end
            end    
            for j=1:size(unique_takeoffs{i},1)
                FN = flt_idxs(find(clust_idxs == clusts(j)));
                for k=1:size(FN,1)
                    plot3(flight_starts_xyz_pruned(1,FN(k)),flight_starts_xyz_pruned(2,FN(k)),flight_starts_xyz_pruned(3,FN(k)),'o','Color',colormap(clusts(j),:),'MarkerSize',8);
                end
            end
        end
    end
    
    % Merge or split takeoff locations if appropriate
    user_input = char(input("If you want to merge clusters, type the two clusters to merge in ascending order (4 8) and the dominant cluster (4). --> '4 8 4'\n"));
    disp('Options are: [] ; "merge X Y X" ; "recluster" \n');
    if isempty(user_input)
        UT = [];
        for i=1:size(unique_takeoffs,2)
            if ~isempty(unique_takeoffs{i})
                UT = [UT,unique_takeoffs{i}+1];
            else
                UT = [UT,NaN];
            end
        end
        UT(1) = 1;
        takeoff_locations_cpu = UT;
        pruned_dataset = non_outlier_indexes;
        % If the clustering is good:
%         save('takeoff_locations_cpu','takeoff_locations_cpu');
%         save('pruned_dataset','pruned_dataset');
%         save('KC','KC');
    elseif strcmp(user_input,'recluster')
        disp("recluster this dataset to get a better split");
        takeoff_locations_cpu = [];
        pruned_dataset = [];
        KC = [];
        
    elseif contains(user_input,"merge") & ~contains(user_input,"split")
        disp("Merging Clusters");
        [KC,unique_takeoffs] = ImBat_MCS_Merge_takeoff_clusters(KC,unique_takeoffs,user_input);   
        UT = [];
        for i=1:size(unique_takeoffs,2)
            if ~isempty(unique_takeoffs{i})
                UT = [UT,unique_takeoffs{i}+1];
            else
                UT = [UT,NaN];
            end
        end
        UT(1) = 1;
        takeoff_locations_cpu = UT;
        pruned_dataset = non_outlier_indexes;

    elseif contains(user_input,"merge") & contains(user_input,"split")
        disp("Splitting & Merging Clusters");
        [KC,unique_takeoffs] = ImBat_MCS_Merge_takeoff_clusters(KC,unique_takeoffs,user_input);   
        UT = [];
        for i=1:size(unique_takeoffs,2)
            if ~isempty(unique_takeoffs{i})
                UT = [UT,unique_takeoffs{i}+1];
            else
                UT = [UT,NaN];
            end
        end
        UT(1) = 1;
        takeoff_locations_cpu = UT;
        pruned_dataset = non_outlier_indexes;        
    elseif ~contains(user_input,"merge") & contains(user_input,"split")
        disp("Splitting Clusters");
        [KC,unique_takeoffs] = ImBat_MCS_Merge_takeoff_clusters(KC,unique_takeoffs,user_input);   
        UT = [];
        for i=1:size(unique_takeoffs,2)
            if ~isempty(unique_takeoffs{i})
                UT = [UT,unique_takeoffs{i}+1];
            else
                UT = [UT,NaN];
            end
        end
        UT(1) = 1;
        takeoff_locations_cpu = UT;
        pruned_dataset = non_outlier_indexes;        
    else
        disp("U didn't type merge or split");
    end
    
    outliers = bin_outlier_flight_indexes;
        
end