% Determine stationarity of transition matrix of flights for different sections of data across training
% 1. Calculate T matrix from the data
% 2. Remove NaN rows
% 3. Simulate a Markov Process
% 4. Recalculate T matrix from the simulated data, check same
% 5. Remove the 1 values and account for them to prevent false absorbing states
% 6. Recalculate T matrix and check for false absorbing states
% 7. Perform test of stationarity 

% Sort the flight starts in time. 
[ss,rr] = sort(flightPaths34.flight_starts_idx);
day_idx = flightPaths34.day(rr);

% Number of different flight types 
li = length(unique(c_s_34));

% Define window size for T matrix calculation
window = 3;
% T is the stack of sliding-window-sized transition matrixes. There will be many of these
% matrixes. 
% The associated flightpaths_in_T_row cell array includes the
% flightpaths that are in each T stack.
T = {}; T_sim = {}; flightpaths_in_T_row = {}; flightpaths_in_T_sim_row = {};
TC_norm_pruned_all = {};
for i=1:size(flightPaths34.Dates,2)-window
    T_temp = zeros(li,li); % Row is F(t), column is F(t+1)
    day_idx_start = find(day_idx==i);
    day_idx_end = find(day_idx==(i+window));
    cs34_range = c_s_34(day_idx_start(1):day_idx_end(end));
    % 1. Calculate T matrix from the data
    for j=1:li
        for m=1:length(cs34_range)-1
            for k=1:li
                if cs34_range(m) == j & cs34_range(m+1) == k
                    T_temp(j,k) =  T_temp(j,k) + 1;
                end
            end
        end
    end
    for j=1:size(T_temp,1)
        row = T_temp(j,:);
        row_t = row./sum(row);
        T_temp_norm(j,:) = row_t;
    end
    % 2. Remove NaN rows
    [T_supp_rows, T_supp_columns] = find(isnan(T_temp_norm));
    [T_supp_rows_nonnan, T_supp_columns_nonnan] = find(~isnan(T_temp_norm));
    flightpaths_in_T_row{i} = unique(T_supp_rows_nonnan);
    T_temp_norm_pruned = T_temp_norm;
    T_temp_norm_pruned(unique(T_supp_rows),:) = [];
    T_temp_norm_pruned(:,unique(T_supp_rows)) = [];
    T{i} = T_temp_norm_pruned;
    
    % 3. Simulate & plot a Markov Process
    clear mc X pre_X
    mc = dtmc(T_temp_norm_pruned);
    numsteps = 1000;
    pre_X = simulate(mc,numsteps);
    % Relabel the things in X
    T_mapping = [unique(T_supp_rows_nonnan),[1:length(unique(T_supp_rows_nonnan))]'];
    for j=1:length(pre_X)
        X(j) = T_mapping(pre_X(j),1);
    end
    %figure; hold on; title(strcat("State Graph of Markov Chain. Days"," ",num2str(i)," ","to"," ",num2str(i+window)));
    %graphplot(mc);
    %figure('name',strcat("Days"," ",num2str(i)," ","to"," ",num2str(i+window))); 
    %simplot(mc,pre_X)
    %figure(); hold on; histogram(X); histogram(cs34_range); title(strcat("Simulated (blue) prior flight ~ and Data (red) prior flight ~ . From window size"," ",num2str(window)));

    % 4. Recalculate T matrix from the simulated data, remove NaN rows
    T_sim_temp = zeros(li,li); % Row is F(t), column is F(t+1)
    for j=1:li
        for m=1:length(X)-1
            for k=1:li
                if X(m) == j & X(m+1) == k
                    T_sim_temp(j,k) =  T_sim_temp(j,k) + 1;
                end
            end
        end
    end
    for j=1:size(T_sim_temp,1)
        row = T_sim_temp(j,:);
        row_t = row./sum(row);
        T_sim_temp_norm(j,:) = row_t;
    end
    
    [T_sim_supp_rows, T_sim_supp_columns] = find(isnan(T_sim_temp_norm));
    [T_sim_supp_rows_nonnan, T_sim_supp_columns_nonnan] = find(~isnan(T_sim_temp_norm));
    flightpaths_in_T_sim_row{i} = unique(T_sim_supp_rows_nonnan);
    T_sim_temp_norm_pruned = T_sim_temp_norm;
    T_sim_temp_norm_pruned(unique(T_sim_supp_rows),:) = [];
    T_sim_temp_norm_pruned(:,unique(T_sim_supp_rows)) = [];
    T_sim{i} = T_sim_temp_norm_pruned;
    
    % Compare T_sim_temp_norm and T_temp_norm_pruned
%     figure();
%     xvalues = string(unique(T_supp_rows_nonnan)');
%     yvalues = string(unique(T_supp_rows_nonnan)');
%     h = heatmap(xvalues,yvalues,T_temp_norm_pruned);
%     h.Title = strcat("Transition Probability Matrix Days",num2str(i)," ","to"," ",num2str(i+1));
%     h.XLabel = 't+1';
%     h.YLabel = 't';
%     figure();
%     xvalues = string(unique(T_supp_rows_nonnan)');
%     yvalues = string(unique(T_supp_rows_nonnan)');
%     h = heatmap(xvalues,yvalues,T_sim_temp_norm_pruned);
%     h.Title = strcat("Simulated Transition Probability Matrix Days",num2str(i)," ","to"," ",num2str(i+1));
%     h.XLabel = 't+1';
%     h.YLabel = 't'; 
    
%     % 5. Remove the 1 row and column from the counts matrix and normalize.
%     clear T_no1_temp_pre T_no1_temp T_no1_temp_norm
%     T_no1_temp_pre = zeros(li,li); % Row is F(t), column is F(t+1)
%     % Calculate the T matrix of this stretch of days 
%     for j=1:li                                            
%         for m=1:length(X)-1
%             for k=1:li
%                 if X(m) == j & X(m+1) == k
%                     T_no1_temp_pre(j,k) =  T_no1_temp_pre(j,k) + 1;
%                 end
%             end
%         end
%     end
%     T_no1_temp = T_no1_temp_pre(2:end,2:end);
%     for j=1:size(T_no1_temp,1)
%         row = T_no1_temp(j,:);
%         row_t = row./sum(row);
%         T_no1_temp_norm(j,:) = row_t;
%     end
%     
%     % Prune the T matrix such that we get a stable state transition matrix
%     % with non-zero entries. INDEXING according to flightpaths_in_T_row
%     [T_no1_supp_rows, T_no1_supp_columns] = find(isnan(T_no1_temp_norm));
%     [T_no1_supp_rows_nonnan, T_no1_supp_columns_nonnan] = find(~isnan(T_no1_temp_norm));
%     temp_flightpaths_in_T_no1_row = unique(T_no1_supp_rows_nonnan);
%     temp_flightpaths_in_T_no1_row(1) = 2; 
%     flightpaths_in_T_no1_row{i} = temp_flightpaths_in_T_no1_row;
%     T_no1_temp_norm_pruned = T_no1_temp_norm;
%     T_no1_temp_norm_pruned(unique(T_no1_supp_rows),:) = [];
%     T_no1_temp_norm_pruned(:,unique(T_no1_supp_rows)) = [];
    
%     figure();
%     xvalues = string(temp_flightpaths_in_T_no1_row');
%     yvalues = string(temp_flightpaths_in_T_no1_row');
%     h = heatmap(xvalues,yvalues,T_no1_temp_norm_pruned);
%     h.Title = strcat("Transition Probability Matrix Days",num2str(i)," ","to"," ",num2str(i+1));
%     h.XLabel = 't+1';
%     h.YLabel = 't';
    
    % Make directed graph of simulated data with no ones. 
%     T_no1_mapping = [unique(T_no1_supp_rows_nonnan),[1:length(unique(T_no1_supp_rows_nonnan))]'];
%     T_no1{i} = T_no1_temp_norm_pruned;
%     clear mc X_no1
%     mc = dtmc(T_no1_temp_norm_pruned);
%     numsteps = 1000;
%     X_no1 = simulate(mc,numsteps);
%     figure; hold on; title(strcat("State Graph of non-1 Markov Chain. Days"," ",num2str(i)," ","to"," ",num2str(i+window)));
%     graphplot(mc);
%     figure('name',strcat("Days"," ",num2str(i)," ","to"," ",num2str(i+window))); 
%     simplot(mc,X_no1)
%     figure(); hold on; histogram(X_no1); histogram(cs34_range); title(strcat("Simulated (blue) prior flight ~ and Data (red) prior flight ~ . No cluster 1. From window size"," ",num2str(window)));
%     
    % 6. That didn't work Michael, so I'll try my method
    [TC_norm_pruned] = ImBat_MCS_helper_weight_ones(cs34_range,li);
    TC_norm_pruned_all{i} = TC_norm_pruned;
    
    
    % Check for false absorbing states
    clear mc X_no1
    mc = dtmc(TC_norm_pruned);
    numsteps = 1000;
    X_no1 = simulate(mc,numsteps);
    figure; hold on; title(strcat("State Graph of non-1 Markov Chain. Days"," ",num2str(i)," ","to"," ",num2str(i+window)));
    graphplot(mc);
    figure('name',strcat("Days"," ",num2str(i)," ","to"," ",num2str(i+window))); 
    simplot(mc,X_no1)
    figure(); hold on; histogram(X_no1); histogram(cs34_range); title(strcat("Simulated (blue) prior flight ~ and Data (red) prior flight ~ . No cluster 1. From window size"," ",num2str(window)));
    
    % 7. Perform test of stationarity? Stationarity on 
end

% End result here is a bunch of 1-weighted matrixes and a bunch of
% non-1-weighted matrixes.




%% WITHOUT CLUSTER 1. SCRAP CODE 
% T is the stack of sliding-window-sized transition matrixes. There will be many of these
% matrixes. The associated flightpaths_in_T_row cell array includes the
% flightpaths that are in each T stack.
T = {}; flightpaths_in_T_row = {};
li = length(unique(c_s_34));
for i=1:size(flightPaths34.Dates,2)-window
    T_temp = zeros(li,li); % Row is F(t), column is F(t+1)
    day_idx_start = find(day_idx==i);
    day_idx_end = find(day_idx==(i+window));
    cs34_range = c_s_34(day_idx_start(1):day_idx_end(end));
    % Calculate the T matrix of this stretch of days 
    for j=2:li                                              % For 2:45
        for m=1:length(cs34_range)-1
            for k=2:li
                if cs34_range(m) == j & cs34_range(m+1) == k
                    T_temp(j,k) =  T_temp(j,k) + 1;
                end
            end
        end
    end
    for j=1:size(T_temp,1)
        row = T_temp(j,:);
        row_t = row./sum(row);
        T_temp_norm(j,:) = row_t;
    end
    % Prune the T matrix such that we get a stable state transition matrix
    % with non-zero entries. INDEXING according to flightpaths_in_T_row
    [T_supp_rows, T_supp_columns] = find(isnan(T_temp_norm));
    [T_supp_rows_nonnan, T_supp_columns_nonnan] = find(~isnan(T_temp_norm));
    flightpaths_in_T_row{i} = unique(T_supp_rows_nonnan);
    T_temp_norm_pruned = T_temp_norm;
    T_temp_norm_pruned(unique(T_supp_rows),:) = [];
    T_temp_norm_pruned(:,unique(T_supp_rows)) = [];
    
    figure();
    xvalues = string(unique(T_supp_rows_nonnan)');
    yvalues = string(unique(T_supp_rows_nonnan)');
    h = heatmap(xvalues,yvalues,T_temp_norm_pruned);
    h.Title = strcat("Transition Probability Matrix Days",num2str(i)," ","to"," ",num2str(i+1));
    h.XLabel = 't+1';
    h.YLabel = 't';
    
    T_mapping = [unique(T_supp_rows_nonnan),[1:length(unique(T_supp_rows_nonnan))]'];
    T{i} = T_temp_norm_pruned;
    clear mc X
    mc = dtmc(T_temp_norm_pruned);
    numsteps = 1000;
    X = simulate(mc,numsteps);
    figure; hold on; title(strcat("State Graph of Markov Chain. Days"," ",num2str(i)," ","to"," ",num2str(i+window)));
    graphplot(mc);
    figure('name',strcat("Days"," ",num2str(i)," ","to"," ",num2str(i+window))); 
    simplot(mc,X)
    figure(); hold on; histogram(X); histogram(cs34_range); title(strcat("Simulated (blue) prior flight ~ and Data (red) prior flight ~ . From window size"," ",num2str(window)));
end



%% DO this for a single day 
% T is the stack of single-day transition matrixes. The associated flightpaths_in_T_row cell array includes the
% flightpaths that are in each T stack.
TT = {}; flightpaths_in_TT_row = {};
li = length(unique(c_s_34));
for i=1:size(flightPaths34.Dates,2)
    T_temp = zeros(li,li); % Row is F(t), column is F(t+1)
    day_idx_start = find(day_idx==i);
    cs34_range = c_s_34(day_idx_start(1):day_idx_start(end));
    % Calculate the T matrix of this stretch of days 
    for j=1:li
        for m=1:length(cs34_range)-1
            for k=1:li
                if cs34_range(m) == j & cs34_range(m+1) == k
                    T_temp(j,k) =  T_temp(j,k) + 1;
                end
            end
        end
    end
    for j=1:size(T_temp,1)
        row = T_temp(j,:);
        row_t = row./sum(row);
        T_temp_norm(j,:) = row_t;
    end
    % Prune the T matrix such that we get a stable state transition matrix
    % with non-zero entries.
    [TT_supp_rows, TT_supp_columns] = find(isnan(T_temp_norm));
    flightpaths_in_TT_row{i} = unique(T_supp_rows);
    T_temp_norm_pruned = T_temp_norm;
    T_temp_norm_pruned(unique(T_supp_rows),:) = [];
    T_temp_norm_pruned(:,unique(T_supp_rows)) = [];
    
    TT{i} = T_temp_norm_pruned;
    clear mc X
    mc = dtmc(T_temp_norm_pruned);
    numsteps = 1000;
    X = simulate(mc,numsteps);
    figure; hold on; title(strcat("State Graph of Markov Chain. Days"," ",num2str(i)," ","to"," ",num2str(i+window)));
    graphplot(mc);
    figure('name',strcat("Days"," ",num2str(i)," ","to"," ",num2str(i+window))); 
    simplot(mc,X)
    figure(); hold on; histogram(X); histogram(cs34_range); title(strcat("Simulated (blue) prior flight ~ and Data (red) prior flight ~ . From window size"," ",num2str(window)));
end













%% For a given window size, calculate the T matrix for moving window
ctr = 0;
KLD_mean_all_windows = {};   
for xx = 4:12%2:floor(length(ROI_Data)/2)-2
    ctr = ctr+1;
    clear T
    li = length(unique(c_s_34));
    window = xx;
    for i=1:size(flightPaths34.Dates,2)-window
        T_temp = zeros(li,li); % Row is F(t), column is F(t+1)
        day_idx_start = find(day_idx==i);
        day_idx_end = find(day_idx==(i+window));
        cs34_range = c_s_34(day_idx_start(1):day_idx_end(end));
        % Calculate the T matrix of this stretch of days 
        for j=1:li
            for m=1:length(cs34_range)-1
                for k=1:li
                    if cs34_range(m) == j & cs34_range(m+1) == k
                        T_temp(j,k) =  T_temp(j,k) + 1;
                    end
                end
            end
        end
        for j=1:size(T_temp,1)
            row = T_temp(j,:);
            row_t = row./sum(row);
            T_temp_norm(j,:) = row_t;
        end
        % Set NaN values to zero 
        T_temp_norm(isnan(T_temp_norm))=0;
        T(:,:,i) = T_temp_norm;          
    end

    % Calculate the mean square error between probability values
    clear MSE MSE_mean
    for i=1:size(T,3)-1
        mse_temp = sqrt((T(:,:,i+1)-T(:,:,i)).^2);
        MSE(i,:) = reshape(mse_temp,[1,size(mse_temp,1)*size(mse_temp,2)])';
        MSE_mean(i) = mean(MSE(i,:));
    end

    %figure(); hold on; plot(MSE_mean); xlabel('Day'); ylabel("MSE between neighboring T matrixes"); title("Stationarity metric");

    %% Calculate the KL divergence between rows 
    clear KLD KLD_mean 
    for i=1:size(T,3)-1
        clear KLD
        Temp_1 = T(:,:,i);
        Temp_2 = T(:,:,i+1);
        for j=1:length(Temp_1)
            Temp_row_1 = Temp_1(j,:);
            Temp_row_2 = Temp_2(j,:);
            if sum(Temp_row_1) ~=1 | sum(Temp_row_2) ~=1
                continue;
            else
                KLD(j) = getKullbackLeibler(Temp_row_1,Temp_row_2);
            end
            if isinf(KLD(j))
                KLD(j) = 1;
            end
        end
        KLD_mean(i) = mean(KLD);
    end
    figure(); subplot(2,1,1); hold on; plot(KLD_mean); plot(MSE_mean*100); xlabel('Day'); ylabel("KLD between neighboring T matrixes"); title(strcat("Stationarity metric Window size "," ",num2str(window)));
    subplot(2,1,2); plot(sig_flight_types);
    KLD_mean_all_windows{ctr} = KLD_mean;
end

for i=1:length(KLD_mean_all_windows)
    [temp,temp_ss] = sort(KLD_mean_all_windows{i});
    min_KL(i) = find(KLD_mean_all_windows{i} == temp(2));
end

% Plot the MSE and the KLD mean 
figure(); hold on; 
for i=1:length(KLD_mean_all_windows); 
    plot(KLD_mean_all_windows{i}); 
end
title(strcat('K-L Divergence of T Matrixes of Sliding Time window '," ",num2str(window))); %plot(MSE_mean);

clear temp;
figure(); hold on;
for i=1:length(KLD_mean_all_windows); 
    temp(i) = mean(KLD_mean_all_windows{i}); 
end
bar(temp); 

% go to ImBat_MCS_twoback_flightplot_PRUNED.m



