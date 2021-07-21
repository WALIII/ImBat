
% Load in the dataset you want 
day_select = 0;
[c_s_34,flightPaths34,num_clust,fd] = ImBat_MCS_Load_Flight_Sequence_Data(ROI_Data,master_track_file,cell_registered_struct,aligned_data_struct,day_select);
Flight_Sequence_Struct.c_s_34 = c_s_34;
Flight_Sequence_Struct.flightpaths = flightPaths34;
Flight_Sequence_Struct.num_clust = num_clust;
Flight_Sequence_Struct.fd = fd;

% Calculate the Transition probability matrix 
[Fnorm34,c_s_Tnorm_34,c_s_T_34,OG_Fnorm] = ImBat_MCS_Calculate_T_Matrix(c_s_34,num_clust,day_select);


%% Create vector of all flight cluster IDs 
true_times = flightPaths34.AllFlightsMasterTime(flightPaths34.flight_starts_idx(:));
[tts,tts_idx]  = sort(true_times(:));
%[r_tts,r_tts_idx] = sort(reward_times(:));
dayx_flight_timeline = flightPaths34.flight_starts_idx(tts_idx)';

% Re-order the flightpaths
positions = sort(flightPaths34.pos(:,:,tts_idx));
total_count = 0;
days_to_examine = max(fd);
prev_day_time = 0;
markov_pst = 0;
R_c_s_34 = [];
for j=1:max(fd)
    day1_indexs = find(fd==j);
    day1_flight_timeline = dayx_flight_timeline(find(fd==j));
    day_1_flights = c_s_34(find(fd==j));
    day1_size = size(day_1_flights,1);
    prev_day_time = day1_flight_timeline(1);
    day1_reward_timeline = ROI_Data{j}.Alignment.out.RewardTime + prev_day_time;
   
    % Timeline Colored Clustered Flights. 
    VCI = [1:max(unique(c_s_34))];
    ll = round([linspace(1,256,size(VCI,2))]);
    colors = gray(size(VCI,2)+12);%hsv(size(VCI,2)+10);
    colorhsvpre = colors(4:end-8,:);
    colorhsv = vertcat([0 0 0],colorhsvpre);
    
    f = figure('Name',strcat('Edited Colored Trajectories and Rewards Day: ',num2str(j)),'NumberTitle','off');
    hold on;
    title(['Timeline of Trajectories']);
    for i=1:size(day1_indexs,1)
        for k=1:size(VCI,2)
            if c_s_34(i+total_count) == VCI(k) 
                xline(day1_flight_timeline(i),'Color',[colorhsv(k,:)],'LineWidth',5);
            end
        end
    end 
    hold off;
    
    f = figure('Name',strcat('3d Trajectories Day: ',num2str(j)),'NumberTitle','off');
    hold on;
    u_flights = unique(c_s_34(1+total_count:size(day1_indexs,1)+total_count));
    for i=1:size(u_flights,1)
        subplot(1,size(u_flights,1),i);
        title(strcat("Cluster ",num2str(u_flights(i))));
        hold on;
        for m=1+total_count:size(day1_indexs,1)+total_count
            if c_s_34(m) == u_flights(i)
                plot3(flightPaths34.pos(1,:,tts_idx(m)),flightPaths34.pos(2,:,tts_idx(m)),flightPaths34.pos(3,:,tts_idx(m)),'Color',[colorhsv(u_flights(i),:)]);
            end
        end
        hold off;
    end  
    
    f = figure('Name',strcat('Raster Timeline Trajectories Day: ',num2str(j)),'NumberTitle','off');
    hold on;
    title(['GuitarHero Plot of Trajectories']);
    ylabel('Flight Cluster');
    for i=1:size(day1_indexs,1)
        for k=1:size(VCI,2)
            if c_s_34(i+total_count) == VCI(k) 
                try
                    line([day1_reward_timeline(i),day1_reward_timeline(i)],[size(VCI,2) 0],'LineStyle','-.','Color',[1.0 0 0]);
                    
                end
                line([day1_flight_timeline(i),day1_flight_timeline(i)],[size(VCI,2) 0],'LineStyle',':','Color',[0.7 0.7 0.7]);
            	line([day1_flight_timeline(i),day1_flight_timeline(i)],[VCI(k)+0.02 VCI(k)+1],'LineWidth',4,'Color',[colorhsv(k,:)]);
            end
        end
    end 
    hold off;
    
    %% Find flights that are rewarded
%     rewarded_flight_indexes = [];
%     for i=1:size(day1_reward_timeline,1)
%         % Find flight that is closest to each reward
%         subtracted_flights = (day1_flight_timeline - day1_reward_timeline(i));
%         closest_flight_index_temp = find(min(subtracted_flights(subtracted_flights>0)));
%         if isempty(closest_flight_index_temp)
%             disp("Reward delivered without preceding flight");
%             rewarded_flight_indexes(i) = 0;
%             continue
%         end
%         closest_flight_idx = closest_flight_index_temp+(size(day1_flight_timeline,1) - size(subtracted_flights(subtracted_flights>0),1)-1);
%         rewarded_flight_indexes(i) = closest_flight_idx;
%     end
%     temp_r = zeros(size(day1_flight_timeline,1),1);
%     nr_idx =nonzeros(rewarded_flight_indexes);
%     temp_r(nr_idx)=1;
%     R_c_s_34 = [R_c_s_34;temp_r];
%     flightPaths34.RewardIDs = R_c_s_34;

% Get inter-flight intervals

        
    if markov_pst == 1
        clear day_markov;
        [out_markov] = ImBat_New_Markov(flightPaths34);
        day_markov.T = out_markov.T;
        day_markov.VA = out_markov.VA(1+total_count:size(day1_indexs,1)+total_count);
        day_markov.out_sort = out_markov.out_sort(1+total_count:size(day1_indexs,1)+total_count);
        day_markov.FlightIDVector = out_markov.FlightIDVector(1+total_count:size(day1_indexs,1)+total_count);
        day_markov.FlightTimeVector = out_markov.FlightTimeVector(1+total_count:size(day1_indexs,1)+total_count);

        % Create prob suff tree
        out_markov.rewardIDs = R_c_s_34;
        ImBat_ProbSuffixTree(out_markov,5); 
    end
    
    total_count = total_count + size(day1_indexs,1);
end
