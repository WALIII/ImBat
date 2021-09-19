% Plot each flight cluster

bat = 'ge';
Cluster = 7;
if bat=='ge'
    sig_struct = [3,4,5,7,12,13,16,17,18];
    no_sig_struct = [1,2,6,8,9,10,11,14,15,19,20];
elseif bat=='ga'
    sig_struct = [5,6,7,8,9];
    no_sig_struct = [1,2,3,4,10,11,12,13,14,15,16];
end

%% Plot each flight type 
flight_height = [0.1, 0.3;
                0.3, 0.5;
                0.5, 0.7;
                0.7, 0.9;
                0.9, 1.1;
                1.1, 1.3;
                1.3, 1.5;
                1.5, 1.7;
                1.7, 1.9;
                1.9, 2.1;
                2.1, 2.3;
                2.3, 2.5;
                2.5, 2.7;
                2.7, 2.9;
                2.9, 3.1;
                3.1, 3.3;
                3.3, 3.5;
                3.5, 3.7;
                3.7, 3.9;
                3.9, 4.1;
                4.1, 4.3;
                4.3, 4.5];

figure(); hold on; title(strcat("All Days; Flight Cluster"," ",num2str(Cluster)));
for r = 1:max(flightPaths34.day)
    day_idx = [];
    for i=1:length(flightPaths34.day)
        if flightPaths34.day(i) == r
            day_idx = [day_idx;i];
        end
    end
    day_flight_starts = flightPaths34.flight_starts_idx(day_idx);
    [ss,rr] = sort(day_flight_starts);
    day_flight_ids = flightPaths34.id(day_idx);
    day_flight_ids_sorted = day_flight_ids(rr);
    day_flight_trajs = flightPaths34.pos(:,:,day_idx);
    day_flight_trajs_sorted = day_flight_trajs(:,:,rr);

    for i=1:length(day_flight_ids_sorted)
        flt_type = day_flight_ids_sorted(i);
        if flt_type == Cluster
            color = [colormap(flt_type,:)];
            plot3(day_flight_trajs_sorted(1,:,i),day_flight_trajs_sorted(2,:,i),day_flight_trajs_sorted(3,:,i),'Color',color);
        end
    end
end

%% Plot the cluster type in each timeline 

for batch = 1:8
    figure(); hold on; title(strcat("Batch"," ",num2str(batch)," ","; All Days; Timeline of only Flight Cluster"," ",num2str(Cluster)));
    for r = (batch-1)*10+1:((batch-1)*10+10) %max(flightPaths34.day)
        day_idx = [];
        for i=1:length(flightPaths34.day)
            if flightPaths34.day(i) == r
                day_idx = [day_idx;i];
            end
        end
        day_flight_starts = flightPaths34.flight_starts_idx(day_idx);
        [ss,rr] = sort(day_flight_starts);
        day_flight_ids = flightPaths34.id(day_idx);
        day_flight_ids_sorted = day_flight_ids(rr);
        day_flight_trajs = flightPaths34.pos(:,:,day_idx);
        day_flight_trajs_sorted = day_flight_trajs(:,:,rr);

        subplot(10,1,r-((batch-1)*10)); hold on;
        for i=1:length(day_flight_ids_sorted)
            flt_type = day_flight_ids_sorted(i);
            if flt_type ~= Cluster 
                color = [0.7 0.7 0.7];
                line([ss(i) ss(i)], [0 5], 'Color',color,'LineStyle',':');
            else
                color = [colormap(flt_type,:)];
                line([ss(i) ss(i)], [0 5], 'Color',color,'LineWidth',2);
            end
        end
    end
end

%% Plot each day's timeline

day=58;
day_idx = [];
for i=1:length(flightPaths34.day)
    if flightPaths34.day(i) == day
        day_idx = [day_idx;i];
    end
end
day_flight_starts = flightPaths34.flight_starts_idx(day_idx);
[ss,rr] = sort(day_flight_starts);
day_flight_ids = flightPaths34.id(day_idx);
day_flight_ids_sorted = day_flight_ids(rr);
day_flight_trajs = flightPaths34.pos(:,:,day_idx);
day_flight_trajs_sorted = day_flight_trajs(:,:,rr);

figure(); hold on; title(strcat("Flight Timeline Day"," ",num2str(day)));
for i=1:length(day_flight_ids_sorted)
    flt_type = day_flight_ids_sorted(i);
%     if ismember(flt_type,sig_struct)
%         scatter(ss(i),flight_height(flt_type,1),'*r');
%     end
    if flt_type >= 20 
        color = [0.7 0.7 0.7];
        line([ss(i) ss(i)], [0 5], 'Color',color,'LineStyle',':');
    elseif ismember(flt_type,sig_struct)
        color = [colormap(flt_type,:)];
        line_height = [flight_height(flt_type,1) flight_height(flt_type,2)];
        line([ss(i) ss(i)], line_height, 'Color',color,'LineWidth',2,'LineStyle',':');
    else
        color = [colormap(flt_type,:)];
        line_height = [flight_height(flt_type,1) flight_height(flt_type,2)];
        line([ss(i) ss(i)], line_height, 'Color',color,'LineWidth',2);
    end
end
% Suspision-- the structure is coming from between bouts     
        
%% Plot each day's timeline in space 

day_idx = [];
for i=1:length(flightPaths34.day)
    if flightPaths34.day(i) == day
        day_idx = [day_idx;i];
    end
end
day_flight_starts = flightPaths34.flight_starts_idx(day_idx);
[ss,rr] = sort(day_flight_starts);
day_flight_ids = flightPaths34.id(day_idx);
day_flight_ids_sorted = day_flight_ids(rr);
day_flight_trajs = flightPaths34.pos(:,:,day_idx);
day_flight_trajs_sorted = day_flight_trajs(:,:,rr);

figure(); hold on; title(strcat("Flight Timeline Day"," ",num2str(day)));
for i=1:length(day_flight_ids_sorted)
    flt_type = day_flight_ids_sorted(i);
    if flt_type >= 20 
        color = [0.7 0.7 0.7];
        plot3(day_flight_trajs_sorted(1,:,i),day_flight_trajs_sorted(2,:,i),day_flight_trajs_sorted(3,:,i),'Color',color,'LineStyle',':');
    else
        color = [colormap(flt_type,:)];
        plot3(day_flight_trajs_sorted(1,:,i),day_flight_trajs_sorted(2,:,i),day_flight_trajs_sorted(3,:,i),'Color',color);
    end
end

%% Does the number of structure flights increase as a function of training?
sig_flight_types = zeros(max(flightPaths34.day),1);
for r = 1:max(flightPaths34.day)
    day_idx = [];
    for i=1:length(flightPaths34.day)
        if flightPaths34.day(i) == r
            day_idx = [day_idx;i];
        end
    end
    day_flight_starts = flightPaths34.flight_starts_idx(day_idx);
    [ss,rr] = sort(day_flight_starts);
    day_flight_ids = flightPaths34.id(day_idx);
    day_flight_ids_sorted = day_flight_ids(rr);
    day_flight_trajs = flightPaths34.pos(:,:,day_idx);
    day_flight_trajs_sorted = day_flight_trajs(:,:,rr);

    for i=1:length(day_flight_ids_sorted)
        flt_type = day_flight_ids_sorted(i);
        if ismember(flt_type,sig_struct)
            sig_flight_types(r) = sig_flight_types(r) + 1;
        end
    end
    sig_flight_types(r) = sig_flight_types(r)/length(day_flight_ids_sorted);
end

[structured_days_most_to_least,order] = sort(sig_flight_types,'descend');
figure(); hold on; plot(sig_flight_types); xlabel("Day"); ylabel("Proportion of P1 Struct FlightPaths"); title("Norm. Sig. Flight Paths across all days");
    

%% Subset the data
%subset_c_s_34 = [];
%subset1_Ge = [14:31];
%subset2_ge = [47:55];
segment = [56:66];
%subset4_ge = [72:76];
subset5_ga = [5,8,12,14];
subset6_ga = [5,7,8,9,12,14,16,18,19,21,23,24,25];
subset_c_s_34 = [];
for r = 1:length(segment)
    day_idx = [];
    for i=1:length(flightPaths34.day)
        if flightPaths34.day(i) == subset6_ga(r)
            day_idx = [day_idx;i];
        end
    end
    day_flight_starts = flightPaths34.flight_starts_idx(day_idx);
    [ss,rr] = sort(day_flight_starts);
    day_flight_ids = flightPaths34.id(day_idx);
    day_flight_ids_sorted = day_flight_ids(rr);
    day_flight_trajs = flightPaths34.pos(:,:,day_idx);
    day_flight_trajs_sorted = day_flight_trajs(:,:,rr);
    subset_c_s_34 = [subset_c_s_34; day_flight_ids_sorted];
end

%% For all days, visualize the flight paths in space to see if you can SEE the stationarity 
for r = 1:max(flightPaths34.day)
    day_idx = [];
    for i=1:length(flightPaths34.day)
        if flightPaths34.day(i) == r
            day_idx = [day_idx;i];
        end
    end
    day_flight_starts = flightPaths34.flight_starts_idx(day_idx);
    [ss,rr] = sort(day_flight_starts);
    day_flight_ids = flightPaths34.id(day_idx);
    day_flight_ids_sorted = day_flight_ids(rr);
    day_flight_trajs = flightPaths34.pos(:,:,day_idx);
    day_flight_trajs_sorted = day_flight_trajs(:,:,rr);

    figure(); hold on; title(strcat("Flight Timeline Day"," ",num2str(r)));
    for i=1:length(day_flight_ids_sorted)
        flt_type = day_flight_ids_sorted(i);
        if flt_type >= 20 
            color = [0.7 0.7 0.7];
            plot3(day_flight_trajs_sorted(1,:,i),day_flight_trajs_sorted(2,:,i),day_flight_trajs_sorted(3,:,i),'Color',color,'LineStyle',':');
        else
            color = [colormap(flt_type,:)];
            plot3(day_flight_trajs_sorted(1,:,i),day_flight_trajs_sorted(2,:,i),day_flight_trajs_sorted(3,:,i),'Color',color);
        end
    end
end

