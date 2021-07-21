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

% 2. Select which subset of days from that dataset you want to include in
% the subsequent analysis (i.e. if there was one bad day, exclude it).
%   day_select is defined in ImBat_MCS_Load_Flight_Sequence_Data.m
%   c_s_34 is the vector of flight types across the entire dataset
%   flightPaths34 is a struct containing flight position & reward info
%   num_clust is the number of unique clusters identified.
%   fd is the vector of what day a given flight belongs to. Same length as c_s_34
day_select = 0;
[c_s_34,flightPaths34,num_clust,fd] = ImBat_MCS_Load_Flight_Sequence_Data(ROI_Data,master_track_file,cell_registered_struct,aligned_data_struct,day_select);

% Calculate the Transition probability matrix (T matrix) 
[Fnorm34,c_s_Tnorm_34,c_s_T_34,OG_Fnorm] = ImBat_MCS_Calculate_T_Matrix(c_s_34,num_clust,day_select);

% Visualize the first-order Markov statistics
sim = 0;
ImBat_MCS_Visualize_DTMC(Fnorm34,c_s_Tnorm_34,c_s_T_34,flightPaths34,sim);

% Simulate sequences from either uniform or weighted (always use weighted
% to make sure the prior flight distribution is like the real data)
weighted = 1;
[legal_sequences_day34,weighted_legal_sequences_day34,WLS_day34,LS_day34,S] = ImBat_MCS_simulated_sequences(Fnorm34,c_s_34,weighted);

% Calculate the Transition probability matrix 
[FnormM,FnormD] = ImBat_MCS_Calculate_Simulated_T_Matrix(S,num_clust,day_select);

% Visualize the distribution of flight types for the data
[flighttype_weights] = ImBat_MCS_Visualize_Prior_Flight_Distribution(c_s_34);

% Visualize the first-order Markov statistics for SIMULATED data
sim=1;
ImBat_MCS_Visualize_DTMC(FnormM,c_s_Tnorm_34,c_s_T_34,flightPaths34,sim);

% Compare the probability of the next flight type conditioned on takeoff position with
% the probability of the next flight type conditioned on takeoff position AND previous flight
% Is there First-Order Structure?
% This function necessittates solving the little Sudoku puzzle of which
% flights takeoff from where. Build a function to do this automatically.
bat_dataset = 'Ge_Longhaul';
[takeoff_locations_cpu] = ImBat_MCS_determine_takeoff_locations(flightPaths34,c_s_34,num_clust,bat_dataset)
ImBat_MCS_Physically_Possible_v_Actual(c_s_34,num_clust)


% [DEPRECATED] Calculate the joint probability of each flight transition
[Joint_OG_Fnorm_vec,Joint_sim_Fnorm_vec] = ImBat_MCS_Joint_probability(OG_Fnorm,FnormM,c_s_34,S);

% [DEPRECATED] Perform a KS test to see if the first order joint probability distributions of the
% real data and simulated data are the same
[h,p,h_row1,p_row1] = ImBat_MCS_KSTest(Joint_sim_Fnorm_vec,Joint_OG_Fnorm_vec,FnormM,OG_Fnorm);


%% Second Order Statistics 

weighted = 1;
[legal_sequences_day34,weighted_legal_sequences_day34,WLS_day34,LS_day34,S] = ImBat_MCS_simulated_sequences(Fnorm34,c_s_34,weighted);

% Calculate the block entropy for 3-grams of the real and simulated data
ngram = 3;
[Entropy,Null_Entropy,Null_Max_Entropy] = ImBat_Block_Entropy(Fnorm34,WLS_day34,LS_day34,weighted_legal_sequences_day34,legal_sequences_day34,c_s_34,day_select,num_clust,ngram);

% Make Jeff's plot
for i=3:6
    ngram = i;
    [Entropy,Null_Entropy,Null_Max_Entropy] = ImBat_Block_Entropy(Fnorm34,WLS_day34,LS_day34,weighted_legal_sequences_day34,legal_sequences_day34,c_s_34,day_select,num_clust,ngram);

    Entropy_ngram(i) = Entropy;
    Null_Entropy_ngram(i) = mean(Null_Entropy);
    Null_Max_Entropy_ngram(i) = mean(Null_Max_Entropy);
end

figure(); hold on; title('Block Entropy'); xlabel('n-gram'); ylabel('Normalized Block Entropy'); plot(Entropy_ngram,'LineWidth',2); plot(Null_Entropy_ngram,'LineWidth',2); plot(Null_Max_Entropy_ngram,'LineWidth',2); 
