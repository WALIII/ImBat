% Workflow for Hippocampus Sequence Project
% MCS 6/7/21

% 0. Put raw data into folders in correct format & Preprocess the raw data
% https://docs.google.com/document/d/1_jYwkI2XeG-4RJJal47YoMkvuga8QrSRvTzrURwCClE/edit)
% This document contains all the info you need to go from raw data to
% dataset with the 4 files listed below.

% 1. Manually load in these 4 files from the prepped dataset
%   aligned_data_struct.mat
%   cellRegistered_....mat
%   Master_Tracking_File_....mat
%   ROI_Data.mat

%   c_s_34 is the vector of flight types across the entire dataset
%   flightPaths34 is a struct containing flight position & reward info
%   num_clust is the number of unique clusters identified.
%   fd is the vector of what day a given flight belongs to. Same length as c_s_34
[c_s_34,flightPaths34,num_clust,fd] = ImBat_MCS_Load_Flight_Sequence_Data(ROI_Data,master_track_file,cell_registered_struct,aligned_data_struct);

% Determine which flight paths have which takeoff locations.
k=6;
[takeoff_locations_cpu,pruned_dataset,KC,outliers] = ImBat_MCS_determine_takeoff_locations(flightPaths34,c_s_34,k);

% Rewind in time to see if the outlier takeoff positions have better
% takeoff positions
rewind_value = 200;
[c_s_r,flightPaths_r,num_clust_r] = ImBat_MCS_rewind_outliers(ROI_Data,master_track_file,cell_registered_struct,aligned_data_struct,outliers,rewind_value);

% Rerun location takeoff determination
k=6;
[takeoff_locations_cpu,pruned_dataset,KC,outliers] = ImBat_MCS_determine_takeoff_locations(flightPaths_r,c_s_r,k);

disp(strcat(num2str(100*(size(outliers,2)/size(c_s_34,1))),"% of the flights of this dataset are outliers"));
% Manually load in KC.mat, takeoff_locations.mat,and pruned_dataset.mat OR
% rerun ImBat_MCS_determine_takeoff_locations before running next step.

% Determine if the conditional probability of the next flight given the
% takeoff location AND the previous flight is different from the
% conditional probability of the next flight given the takeoff location,
% previous flight AND previous takeoff location. THis tells us if the
% actual flight they take from A to B gives us information about what kind
% of flight they will fly next, while holding constant the locations of
% where they come from.
[h,p] = ImBat_MCS_TwoBack_Physically_Possible_v_Actual(c_s_34,pruned_dataset,takeoff_locations_cpu,outliers);

% Calculate the Transition probability matrix (T matrix) 
[Fnorm,c_s_Tnorm,c_s_T,OG_Fnorm] = ImBat_MCS_Calculate_T_Matrix(c_s_34,num_clust);

% Visualize the first-order Markov statistics
sim=0;
ImBat_MCS_Visualize_DTMC(Fnorm,c_s_Tnorm,c_s_T,flightPaths34,sim);

% Simulate sequences from either uniform or weighted 
weighted = 1;
[legal_sequences_day,weighted_legal_sequences_day,WLS_day,LS_day,S] = ImBat_MCS_simulated_sequences(Fnorm,c_s_34,weighted);

% Determine if the conditional probability of the next flight given the takeoff location 
% is different from the conditional probability of the next flight taking into account
% the takeoff location AND the previous flights
sim_outliers = [];
[H_SIM,P_SIM] = ImBat_MCS_TwoBack_Physically_Possible_v_Actual_simulated(S,pruned_dataset,takeoff_locations_cpu,sim_outliers);

% See how many S are equal to c_s_34
isequal_tally = 0;
for i=1:size(S,1)
    if isequal(S(i,:)',c_s_34)
        isequal_tally = isequal_tally+1;
    end
end

% Calculate the Transition probability matrix 
[FnormM,FnormD] = ImBat_MCS_Calculate_Simulated_T_Matrix(S,num_clust,day_select);

% Visualize the first-order Markov statistics for SIMULATED data
sim=1;
ImBat_MCS_Visualize_DTMC(FnormM,c_s_Tnorm,c_s_T,flightPaths34,sim);

% Calculate the joint probability of each flight transition
% [Joint_OG_Fnorm_vec,Joint_sim_Fnorm_vec] = ImBat_MCS_Joint_probability(OG_Fnorm,FnormM,c_s_34,S);
% 
% % Perform a KS test to see if the first order joint probability distributions of the
% % real data and simulated data are the same
% [h,p,h_row1,p_row1] = ImBat_MCS_KSTest(Joint_sim_Fnorm_vec,Joint_OG_Fnorm_vec,FnormM,OG_Fnorm);


%% Second Order Statistics 

% Calculate the block entropy for 3-grams of the real and simulated data
ngram = 3;
[Entropy,Null_Entropy,Null_Max_Entropy] = ImBat_Block_Entropy(Fnorm,WLS_day,LS_day,weighted_legal_sequences_day,legal_sequences_day,c_s_34,day_select,num_clust,ngram);

% Make Jeff's plot
for i=3:6
    ngram = i;
    [Entropy,Null_Entropy,Null_Max_Entropy] = ImBat_Block_Entropy(Fnorm,WLS_day,LS_day,weighted_legal_sequences_day,legal_sequences_day,c_s_34,day_select,num_clust,ngram);

    Entropy_ngram(i) = Entropy;
    Null_Entropy_ngram(i) = mean(Null_Entropy);
    Null_Max_Entropy_ngram(i) = mean(Null_Max_Entropy);
end

figure(); hold on; title('Block Entropy'); xlabel('n-gram'); ylabel('Normalized Block Entropy'); plot(Entropy_ngram,'LineWidth',2); plot(Null_Entropy_ngram,'LineWidth',2); plot(Null_Max_Entropy_ngram,'LineWidth',2); 
