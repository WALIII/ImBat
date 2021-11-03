% Counting how much data I have

minimum_number_of_transition_instances = 20;

%% Genarro 
disp("Beginning Genarro Stats")
ge_c_s_34 = load('Processed_Ge_LongHaul/c_s_34.mat');
ge_flightPaths = load('Processed_Ge_LongHaul/FlightPaths34.mat');
[co,number_unique_transitions_2,number_unique_transitions_3] = ImBat_MCS_Count_Unique_Transitions(ge_c_s_34.c_s_34,ge_flightPaths.flightPaths34,minimum_number_of_transition_instances)

% Log data in grand total struct 
grand_totals.ge.number_unique_transitions_2 = number_unique_transitions_2;
grand_totals.ge.number_unique_transitions_3 = number_unique_transitions_3;
grand_totals.ge.number_flights = length(ge_c_s_34.c_s_34);

clear ge_flightPaths ge_c_s_34;

% Examine subsets of the data to look for how many examples of a given
% transition type we have in a window that a cell can be aligned across
subset_type=1;
if subset_type == 1
load('Processed_Ge_5_thru_11/CellReg_files/ROI_Data/cellRegistered_20210110_211744.mat');
load('Processed_Ge_5_thru_11/CellReg_files/ROI_Data/aligned_data_struct.mat');
load('Processed_Ge_5_thru_11/CellReg_files/ROI_Data/Master_Tracking_File_19-Nov-2020 08:59:12.mat');
load('Processed_Ge_5_thru_11/CellReg_files/ROI_Data/ROI_Data.mat');
elseif subset_type == 2
load('Processed_Ge_19_thru_24/aligned_data_struct.mat');
load('Processed_Ge_19_thru_24/cellRegistered_20210126_202011.mat');
load('Processed_Ge_19_thru_24/Master_Tracking_File_19-Nov-2020 08:59:12.mat');
load('Processed_Ge_19_thru_24/ROI_Data.mat');
end
[CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI] = ImBat_MCS_Align_FlightPaths(aligned_data_struct,cell_registered_struct,master_track_file,ROI_Data,co);

cell_set = 1; %[1, 2, 6, 13, 17, 20, 25, 26, 27, 37, 38, 41, 43, 47, 48, 67, 69, 93];
day_set = 1:5;
TG_cat = []; DG_cat = [];
for i=1:co
    cluster_number = i;
    [TG_transition_number_list,DG_transition_number_list] = ImBat_MCS_Unique_Transition_Splits_Subsetted_Data(CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI,cell_set,day_set,cluster_number);
    TG_transition_agg{i} = TG_transition_number_list;
    DG_transition_agg{i} = DG_transition_number_list;
    TG_cat = [TG_cat;TG_transition_number_list];
    DG_cat = [DG_cat;DG_transition_number_list];
end

% How many transitions with at least 10 examples are there 
TG_greater_than_10 = length(find(TG_cat(TG_cat>=10)));
DG_greater_than_10 = length(find(DG_cat(DG_cat>=10)));

% How many transitions with at least 30 examples are there
TG_greater_than_30 = length(find(TG_cat(TG_cat>=30)));
DG_greater_than_30 = length(find(DG_cat(DG_cat>=30)));

% Log totals
grand_totals.ge.subset1.number_unique_TG_GT10 = TG_greater_than_10;
grand_totals.ge.subset1.number_unique_DG_GT10 = DG_greater_than_10;
grand_totals.ge.subset1.number_unique_TG_GT30 = TG_greater_than_30;
grand_totals.ge.subset1.number_unique_DG_GT30 = DG_greater_than_30;

% Display Stats
disp("Genarro");
disp(strcat("Diversity and Quantity of Flight Behavior over"," ",num2str(length(ge_flightPaths.flightPaths34.Dates))," ","days of training."));
disp(strcat("Total # Flights:"," ",num2str(length(ge_c_s_34.c_s_34))));
disp(strcat("Total # Unique Transition Types (2-gram):"," ",num2str(number_unique_transitions_2)));
disp(strcat("Total # Unique Transition Types (3-gram):"," ",num2str(number_unique_transitions_3)));
disp(strcat("Total # of 3-gram transition types with >=10 occurances:"," ",num2str(TG_greater_than_10)));
disp(strcat("Total # of 2-gram transition types with >=10 occurances:"," ",num2str(DG_greater_than_10)));
disp(strcat("Total # of 3-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));
disp(strcat("Total # of 2-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));

clear ROI_Data aligned_data_struct cell_registered_struct master_track_file;
clear TG_greater_than_10 DG_greater_than_10 DG_greater_than_30 TG_greater_than_30;
clear TG_transition_agg DG_transition_agg TG_transition_number_list DG_transition_number_list;
clear CombinedROI ROI_Data FlightPaths c_s_34 FlightAlignedROI;
clear co number_unique_transitions_2 number_unique_transitions_3;

%% Galileo 
close all;
disp("Beginning Galileo Stats");
ga_c_s_34 = load('Processed_Ga_LongHaul/c_s_34.mat');
ga_flightPaths = load('Processed_Ga_LongHaul/flightPaths34.mat');
[co,number_unique_transitions_2,number_unique_transitions_3] = ImBat_MCS_Count_Unique_Transitions(ga_c_s_34.c_s_34,ga_flightPaths.flightPaths34,minimum_number_of_transition_instances)

% Log data in grand total struct 
grand_totals.ga.number_unique_transitions_2 = number_unique_transitions_2;
grand_totals.ga.number_unique_transitions_3 = number_unique_transitions_3;
grand_totals.ga.number_flights = length(ga_c_s_34.c_s_34);

num_days = length(ga_flightPaths.flightPaths34.Dates);
num_flights = length(ga_c_s_34.c_s_34);
clear ga_flightPaths ga_c_s_34;

% Examine subsets of the data to look for how many examples of a given
% transition type we have in a window that a cell can be aligned across
load('Processed_Ga_11_thru_24/cellRegistered_20210203_111252.mat');
load('Processed_Ga_11_thru_24/aligned_data_struct.mat');
load('Processed_Ga_11_thru_24/Master_Tracking_File_19-Nov-2020 08:59:12.mat');
load('Processed_Ga_11_thru_24/ROI_Data.mat');
[CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI] = ImBat_MCS_Align_FlightPaths(aligned_data_struct,cell_registered_struct,master_track_file,ROI_Data,co);

cell_set = 1; %[1, 2, 6, 13, 17, 20, 25, 26, 27, 37, 38, 41, 43, 47, 48, 67, 69, 93];
day_set = 1:5;
TG_cat = []; DG_cat = [];
for i=1:co
    cluster_number = i;
    [TG_transition_number_list,DG_transition_number_list] = ImBat_MCS_Unique_Transition_Splits_Subsetted_Data(CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI,cell_set,day_set,cluster_number);
    TG_transition_agg{i} = TG_transition_number_list;
    DG_transition_agg{i} = DG_transition_number_list;
    TG_cat = [TG_cat;TG_transition_number_list];
    DG_cat = [DG_cat;DG_transition_number_list];
end

% How many transitions with at least 10 examples are there 
TG_greater_than_10 = length(find(TG_cat(TG_cat>=10)));
DG_greater_than_10 = length(find(DG_cat(DG_cat>=10)));

% How many transitions with at least 30 examples are there
TG_greater_than_30 = length(find(TG_cat(TG_cat>=30)));
DG_greater_than_30 = length(find(DG_cat(DG_cat>=30)));

% Log totals
grand_totals.ga.subset1.number_unique_TG_GT10 = TG_greater_than_10;
grand_totals.ga.subset1.number_unique_DG_GT10 = DG_greater_than_10;
grand_totals.ga.subset1.number_unique_TG_GT30 = TG_greater_than_30;
grand_totals.ga.subset1.number_unique_DG_GT30 = DG_greater_than_30;

% Display Stats
disp("Galileo");
disp(strcat("Diversity and Quantity of Flight Behavior over"," ",num2str(num_days)," ","days of training."));
disp(strcat("Total # Flights:"," ",num2str(num_flights)));
disp(strcat("Total # Unique Transition Types (2-gram):"," ",num2str(number_unique_transitions_2)));
disp(strcat("Total # Unique Transition Types (3-gram):"," ",num2str(number_unique_transitions_3)));
disp(strcat("Total # of 3-gram transition types with >=10 occurances:"," ",num2str(TG_greater_than_10)));
disp(strcat("Total # of 2-gram transition types with >=10 occurances:"," ",num2str(DG_greater_than_10)));
disp(strcat("Total # of 3-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));
disp(strcat("Total # of 2-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));

clear ROI_Data aligned_data_struct cell_registered_struct master_track_file;
clear TG_greater_than_10 DG_greater_than_10 DG_greater_than_30 TG_greater_than_30;
clear TG_transition_agg DG_transition_agg TG_transition_number_list DG_transition_number_list;
clear CombinedROI ROI_Data FlightPaths c_s_34 FlightAlignedROI;
clear co number_unique_transitions_2 number_unique_transitions_3;

%% Gongolo
close all;
disp("Beginning Gongolo Stats");
go_c_s_34 = load('Processed_Go_26_thru_03/CellReg_files/ROI_Data/c_s_34.mat');
go_flightPaths = load('Processed_Go_26_thru_03/CellReg_files/ROI_Data/flightPaths34.mat');
[co,number_unique_transitions_2,number_unique_transitions_3] = ImBat_MCS_Count_Unique_Transitions(go_c_s_34.c_s_34,go_flightPaths.flightPaths34,minimum_number_of_transition_instances)

% Log data in grand total struct 
grand_totals.go.number_unique_transitions_2 = number_unique_transitions_2;
grand_totals.go.number_unique_transitions_3 = number_unique_transitions_3;
grand_totals.go.number_flights = length(go_c_s_34.c_s_34);

num_days = length(go_flightPaths.flightPaths34.Dates);
num_flights = length(go_c_s_34.c_s_34);
clear go_flightPaths go_c_s_34;

% Examine subsets of the data to look for how many examples of a given
% transition type we have in a window that a cell can be aligned across
load('Processed_Go_26_thru_03/CellReg_files/ROI_Data/aligned_data_struct.mat');
load('Processed_Go_26_thru_03/CellReg_files/ROI_Data/cellRegistered_20201104_203847.mat');
load('Processed_Go_26_thru_03/CellReg_files/ROI_Data/Master_Tracking_File.mat');
load('Processed_Go_26_thru_03/CellReg_files/ROI_Data/ROI_Data.mat');
[CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI] = ImBat_MCS_Align_FlightPaths(aligned_data_struct,cell_registered_struct,master_track_file,ROI_Data,co);

cell_set = 1; %[1, 2, 6, 13, 17, 20, 25, 26, 27, 37, 38, 41, 43, 47, 48, 67, 69, 93];
day_set = 1:5;
TG_cat = []; DG_cat = [];
for i=1:co
    cluster_number = i;
    [TG_transition_number_list,DG_transition_number_list] = ImBat_MCS_Unique_Transition_Splits_Subsetted_Data(CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI,cell_set,day_set,cluster_number);
    TG_transition_agg{i} = TG_transition_number_list;
    DG_transition_agg{i} = DG_transition_number_list;
    TG_cat = [TG_cat;TG_transition_number_list];
    DG_cat = [DG_cat;DG_transition_number_list];
end

% How many transitions with at least 10 examples are there 
TG_greater_than_10 = length(find(TG_cat(TG_cat>=10)));
DG_greater_than_10 = length(find(DG_cat(DG_cat>=10)));

% How many transitions with at least 30 examples are there
TG_greater_than_30 = length(find(TG_cat(TG_cat>=30)));
DG_greater_than_30 = length(find(DG_cat(DG_cat>=30)));

% Log totals
grand_totals.go.subset1.number_unique_TG_GT10 = TG_greater_than_10;
grand_totals.go.subset1.number_unique_DG_GT10 = DG_greater_than_10;
grand_totals.go.subset1.number_unique_TG_GT30 = TG_greater_than_30;
grand_totals.go.subset1.number_unique_DG_GT30 = DG_greater_than_30;

% Display Stats
disp("Gongolo");
disp(strcat("Diversity and Quantity of Flight Behavior over"," ",num2str(num_days)," ","days of training."));
disp(strcat("Total # Flights:"," ",num2str(num_flights)));
disp(strcat("Total # Unique Transition Types (2-gram):"," ",num2str(number_unique_transitions_2)));
disp(strcat("Total # Unique Transition Types (3-gram):"," ",num2str(number_unique_transitions_3)));
disp(strcat("Total # of 3-gram transition types with >=10 occurances:"," ",num2str(TG_greater_than_10)));
disp(strcat("Total # of 2-gram transition types with >=10 occurances:"," ",num2str(DG_greater_than_10)));
disp(strcat("Total # of 3-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));
disp(strcat("Total # of 2-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));

clear ROI_Data aligned_data_struct cell_registered_struct master_track_file;
clear TG_greater_than_10 DG_greater_than_10 DG_greater_than_30 TG_greater_than_30;
clear TG_transition_agg DG_transition_agg TG_transition_number_list DG_transition_number_list;
clear CombinedROI ROI_Data FlightPaths c_s_34 FlightAlignedROI;
clear co number_unique_transitions_2 number_unique_transitions_3;

%% Giovanni
disp("Beginning Giovanni Stats");
gi_c_s_34 = load('Processed_Gi_07_thru_18/c_s_34.mat');
gi_flightPaths = load('Processed_Gi_07_thru_18/flightPaths34.mat');
[co,number_unique_transitions_2,number_unique_transitions_3] = ImBat_MCS_Count_Unique_Transitions(gi_c_s_34.c_s_34,gi_flightPaths.flightPaths34,minimum_number_of_transition_instances)

% Log data in grand total struct 
grand_totals.gi.number_unique_transitions_2 = number_unique_transitions_2;
grand_totals.gi.number_unique_transitions_3 = number_unique_transitions_3;
grand_totals.gi.number_flights = length(gi_c_s_34.c_s_34);

num_days = length(gi_flightPaths.flightPaths34.Dates);
num_flights = length(gi_c_s_34.c_s_34);
clear gi_flightPaths gi_c_s_34;

% Examine subsets of the data to look for how many examples of a given
% transition type we have in a window that a cell can be aligned across
load('Processed_Gi_07_thru_18/cellRegistered_20201216_184647.mat');
load('Processed_Gi_07_thru_18/aligned_data_struct.mat');
load('Processed_Gi_07_thru_18/ROI_Data.mat');
load('Processed_Gi_07_thru_18/Master_Tracking_File.mat');
[CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI] = ImBat_MCS_Align_FlightPaths(aligned_data_struct,cell_registered_struct,master_track_file,ROI_Data,co);

cell_set = 1; %[1, 2, 6, 13, 17, 20, 25, 26, 27, 37, 38, 41, 43, 47, 48, 67, 69, 93];
day_set = 1:5;
TG_cat = []; DG_cat = [];
for i=1:co
    cluster_number = i;
    [TG_transition_number_list,DG_transition_number_list] = ImBat_MCS_Unique_Transition_Splits_Subsetted_Data(CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI,cell_set,day_set,cluster_number);
    TG_transition_agg{i} = TG_transition_number_list;
    DG_transition_agg{i} = DG_transition_number_list;
    TG_cat = [TG_cat;TG_transition_number_list];
    DG_cat = [DG_cat;DG_transition_number_list];
end

% How many transitions with at least 10 examples are there 
TG_greater_than_10 = length(find(TG_cat(TG_cat>=10)));
DG_greater_than_10 = length(find(DG_cat(DG_cat>=10)));

% How many transitions with at least 30 examples are there
TG_greater_than_30 = length(find(TG_cat(TG_cat>=30)));
DG_greater_than_30 = length(find(DG_cat(DG_cat>=30)));

% Log totals
grand_totals.gi.subset1.number_unique_TG_GT10 = TG_greater_than_10;
grand_totals.gi.subset1.number_unique_DG_GT10 = DG_greater_than_10;
grand_totals.gi.subset1.number_unique_TG_GT30 = TG_greater_than_30;
grand_totals.gi.subset1.number_unique_DG_GT30 = DG_greater_than_30;

% Display Stats
disp("Giovanni");
disp(strcat("Diversity and Quantity of Flight Behavior over"," ",num2str(num_days)," ","days of training."));
disp(strcat("Total #o Flights:"," ",num2str(num_flights)));
disp(strcat("Total # Unique Transition Types (2-gram):"," ",num2str(number_unique_transitions_2)));
disp(strcat("Total # Unique Transition Types (3-gram):"," ",num2str(number_unique_transitions_3)));
disp(strcat("Total # of 3-gram transition types with >=10 occurances:"," ",num2str(TG_greater_than_10)));
disp(strcat("Total # of 2-gram transition types with >=10 occurances:"," ",num2str(DG_greater_than_10)));
disp(strcat("Total # of 3-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));
disp(strcat("Total # of 2-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));

clear ROI_Data aligned_data_struct cell_registered_struct master_track_file;
clear TG_greater_than_10 DG_greater_than_10 DG_greater_than_30 TG_greater_than_30;
clear TG_transition_agg DG_transition_agg TG_transition_number_list DG_transition_number_list;
clear CombinedROI ROI_Data FlightPaths c_s_34 FlightAlignedROI;
clear co number_unique_transitions_2 number_unique_transitions_3;

%% Zack
close all;
disp("Beginning Zack Stats");
za_c_s_34 = load('Processed_za_LongHaul/CellReg_files/ROI_Data/c_s_34.mat');
za_flightPaths = load('Processed_za_LongHaul/CellReg_files/ROI_Data/flightPaths34.mat');
[co,number_unique_transitions_2,number_unique_transitions_3] = ImBat_MCS_Count_Unique_Transitions(za_c_s_34.c_s_34,za_flightPaths.flightPaths34,minimum_number_of_transition_instances)

% Log data in grand total struct 
grand_totals.za.number_unique_transitions_2 = number_unique_transitions_2;
grand_totals.za.number_unique_transitions_3 = number_unique_transitions_3;
grand_totals.za.number_flights = length(za_c_s_34.c_s_34);

num_flights = length(za_c_s_34.c_s_34);
num_days = length(za_flightPaths.flightPaths34.Dates);
clear za_flightPaths za_c_s_34;

% Examine subsets of the data to look for how many examples of a given
% transition type we have in a window that a cell can be aligned across
load('Processed_za_LongHaul/CellReg_files/ROI_Data/aligned_data_struct.mat');
load('Processed_za_LongHaul/CellReg_files/ROI_Data/cellRegistered_20210107_090526.mat');
load('Processed_za_LongHaul/CellReg_files/ROI_Data/ROI_Data.mat');
load('Processed_za_LongHaul/CellReg_files/ROI_Data/Master_Tracking_File_19-Nov-2020 08:59:12.mat');
[CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI] = ImBat_MCS_Align_FlightPaths(aligned_data_struct,cell_registered_struct,master_track_file,ROI_Data,co);

cell_set = 1; %[1, 2, 6, 13, 17, 20, 25, 26, 27, 37, 38, 41, 43, 47, 48, 67, 69, 93];
day_set = 1:5;
TG_cat = []; DG_cat = [];
for i=1:co
    cluster_number = i;
    [TG_transition_number_list,DG_transition_number_list] = ImBat_MCS_Unique_Transition_Splits_Subsetted_Data(CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI,cell_set,day_set,cluster_number);
    TG_transition_agg{i} = TG_transition_number_list;
    DG_transition_agg{i} = DG_transition_number_list;
    TG_cat = [TG_cat;TG_transition_number_list];
    DG_cat = [DG_cat;DG_transition_number_list];
end

% How many transitions with at least 10 examples are there 
TG_greater_than_10 = length(find(TG_cat(TG_cat>=10)));
DG_greater_than_10 = length(find(DG_cat(DG_cat>=10)));

% How many transitions with at least 30 examples are there
TG_greater_than_30 = length(find(TG_cat(TG_cat>=30)));
DG_greater_than_30 = length(find(DG_cat(DG_cat>=30)));

% Log totals
grand_totals.za.subset1.number_unique_TG_GT10 = TG_greater_than_10;
grand_totals.za.subset1.number_unique_DG_GT10 = DG_greater_than_10;
grand_totals.za.subset1.number_unique_TG_GT30 = TG_greater_than_30;
grand_totals.za.subset1.number_unique_DG_GT30 = DG_greater_than_30;

% Display Stats
disp("Zack");
disp(strcat("Diversity and Quantity of Flight Behavior over"," ",num2str(num_days)," ","days of training."));
disp(strcat("Total # Flights:"," ",num2str(num_flights)));
disp(strcat("Total # Unique Transition Types (2-gram):"," ",num2str(number_unique_transitions_2)));
disp(strcat("Total # Unique Transition Types (3-gram):"," ",num2str(number_unique_transitions_3)));
disp(strcat("Total # of 3-gram transition types with >=10 occurances:"," ",num2str(TG_greater_than_10)));
disp(strcat("Total # of 2-gram transition types with >=10 occurances:"," ",num2str(DG_greater_than_10)));
disp(strcat("Total # of 3-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));
disp(strcat("Total # of 2-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));

clear ROI_Data aligned_data_struct cell_registered_struct master_track_file;
clear TG_greater_than_10 DG_greater_than_10 DG_greater_than_30 TG_greater_than_30;
clear TG_transition_agg DG_transition_agg TG_transition_number_list DG_transition_number_list;
clear CombinedROI ROI_Data FlightPaths c_s_34 FlightAlignedROI;
clear co number_unique_transitions_2 number_unique_transitions_3;

%% Zuzu
close all;
disp("Beginning Zuzu Stats");
zu_c_s_34 = load('Processed_zu_09_thru_20/CellReg_files/ROI_Data/c_s_34.mat');
zu_flightPaths = load('Processed_zu_09_thru_20/CellReg_files/ROI_Data/flightPaths.mat');
[co,number_unique_transitions_2,number_unique_transitions_3] = ImBat_MCS_Count_Unique_Transitions(zu_c_s_34.c_s_34,zu_flightPaths.flightPaths,minimum_number_of_transition_instances)

% Log data in grand total struct 
grand_totals.zu.number_unique_transitions_2 = number_unique_transitions_2;
grand_totals.zu.number_unique_transitions_3 = number_unique_transitions_3;
grand_totals.zu.number_flights = length(zu_c_s_34.c_s_34);

num_flights = length(zu_c_s_34.c_s_34);
num_days = length(zu_flightPaths.flightPaths.Dates);
clear zu_flightPaths zu_c_s_34;

% Examine subsets of the data to look for how many examples of a given
% transition type we have in a window that a cell can be aligned across
load('Processed_zu_09_thru_20/CellReg_files/ROI_Data/aligned_data_struct.mat');
load('Processed_zu_09_thru_20/CellReg_files/ROI_Data/cellRegistered_20201217_131238.mat');
load('Processed_zu_09_thru_20/CellReg_files/ROI_Data/ROI_Data.mat');
load('Processed_zu_09_thru_20/CellReg_files/ROI_Data/Master_Tracking_File.mat');
[CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI] = ImBat_MCS_Align_FlightPaths(aligned_data_struct,cell_registered_struct,master_track_file,ROI_Data,co);

cell_set = 1; %[1, 2, 6, 13, 17, 20, 25, 26, 27, 37, 38, 41, 43, 47, 48, 67, 69, 93];
day_set = 1:5;
TG_cat = []; DG_cat = [];
for i=1:co
    cluster_number = i;
    [TG_transition_number_list,DG_transition_number_list] = ImBat_MCS_Unique_Transition_Splits_Subsetted_Data(CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI,cell_set,day_set,cluster_number);
    TG_transition_agg{i} = TG_transition_number_list;
    DG_transition_agg{i} = DG_transition_number_list;
    TG_cat = [TG_cat;TG_transition_number_list];
    DG_cat = [DG_cat;DG_transition_number_list];
end

% How many transitions with at least 10 examples are there 
TG_greater_than_10 = length(find(TG_cat(TG_cat>=10)));
DG_greater_than_10 = length(find(DG_cat(DG_cat>=10)));

% How many transitions with at least 30 examples are there
TG_greater_than_30 = length(find(TG_cat(TG_cat>=30)));
DG_greater_than_30 = length(find(DG_cat(DG_cat>=30)));

% Log totals
grand_totals.zu.subset1.number_unique_TG_GT10 = TG_greater_than_10;
grand_totals.zu.subset1.number_unique_DG_GT10 = DG_greater_than_10;
grand_totals.zu.subset1.number_unique_TG_GT30 = TG_greater_than_30;
grand_totals.zu.subset1.number_unique_DG_GT30 = DG_greater_than_30;

% Display Stats
disp("Zuzu");
disp(strcat("Diversity and Quantity of Flight Behavior over"," ",num2str(num_days)," ","days of training."));
disp(strcat("Total # Flights:"," ",num2str(num_flights)));
disp(strcat("Total # Unique Transition Types (2-gram):"," ",num2str(number_unique_transitions_2)));
disp(strcat("Total # Unique Transition Types (3-gram):"," ",num2str(number_unique_transitions_3)));
disp(strcat("Total # of 3-gram transition types with >=10 occurances:"," ",num2str(TG_greater_than_10)));
disp(strcat("Total # of 2-gram transition types with >=10 occurances:"," ",num2str(DG_greater_than_10)));
disp(strcat("Total # of 3-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));
disp(strcat("Total # of 2-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));

clear ROI_Data aligned_data_struct cell_registered_struct master_track_file;
clear TG_greater_than_10 DG_greater_than_10 DG_greater_than_30 TG_greater_than_30;
clear TG_transition_agg DG_transition_agg TG_transition_number_list DG_transition_number_list;
clear CombinedROI ROI_Data FlightPaths c_s_34 FlightAlignedROI;
clear co number_unique_transitions_2 number_unique_transitions_3;

%% Zoom2
close all;
disp("Beginning Zoom2 Stats");
z2_c_s_34 = load('Processed_z2_LongHaul/CellReg_files/ROI_Data/c_s_34.mat');
z2_flightPaths = load('Processed_z2_LongHaul/CellReg_files/ROI_Data/flightPaths34.mat');
[co,number_unique_transitions_2,number_unique_transitions_3] = ImBat_MCS_Count_Unique_Transitions(z2_c_s_34.c_s_34,z2_flightPaths.flightPaths34,minimum_number_of_transition_instances)

% Log data in grand total struct 
grand_totals.z2.number_unique_transitions_2 = number_unique_transitions_2;
grand_totals.z2.number_unique_transitions_3 = number_unique_transitions_3;
grand_totals.z2.number_flights = length(z2_c_s_34.c_s_34);

num_flights = length(z2_c_s_34.c_s_34);
num_days = length(z2_flightPaths.flightPaths34.Dates);
clear z2_flightPaths z2_c_s_34;

% Examine subsets of the data to look for how many examples of a given
% transition type we have in a window that a cell can be aligned across
load('Processed_z2_LongHaul/CellReg_files/ROI_Data/aligned_data_struct.mat');
load('Processed_z2_LongHaul/CellReg_files/ROI_Data/cellRegistered_20210112_074352.mat');
load('Processed_z2_LongHaul/CellReg_files/ROI_Data/Master_Tracking_File_19-Nov-2020 08:59:12.mat');
load('Processed_z2_LongHaul/CellReg_files/ROI_Data/ROI_Data.mat');
[CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI] = ImBat_MCS_Align_FlightPaths(aligned_data_struct,cell_registered_struct,master_track_file,ROI_Data,co);

cell_set = 1; %[1, 2, 6, 13, 17, 20, 25, 26, 27, 37, 38, 41, 43, 47, 48, 67, 69, 93];
day_set = 1:5;
TG_cat = []; DG_cat = [];
for i=1:co
    cluster_number = i;
    [TG_transition_number_list,DG_transition_number_list] = ImBat_MCS_Unique_Transition_Splits_Subsetted_Data(CombinedROI,ROI_Data,FlightPaths,c_s_34,FlightAlignedROI,cell_set,day_set,cluster_number);
    TG_transition_agg{i} = TG_transition_number_list;
    DG_transition_agg{i} = DG_transition_number_list;
    TG_cat = [TG_cat;TG_transition_number_list];
    DG_cat = [DG_cat;DG_transition_number_list];
end

% How many transitions with at least 10 examples are there 
TG_greater_than_10 = length(find(TG_cat(TG_cat>=10)));
DG_greater_than_10 = length(find(DG_cat(DG_cat>=10)));

% How many transitions with at least 30 examples are there
TG_greater_than_30 = length(find(TG_cat(TG_cat>=30)));
DG_greater_than_30 = length(find(DG_cat(DG_cat>=30)));

% Log totals
grand_totals.z2.subset1.number_unique_TG_GT10 = TG_greater_than_10;
grand_totals.z2.subset1.number_unique_DG_GT10 = DG_greater_than_10;
grand_totals.z2.subset1.number_unique_TG_GT30 = TG_greater_than_30;
grand_totals.z2.subset1.number_unique_DG_GT30 = DG_greater_than_30;

% Display Stats
disp("Zoom2");
disp(strcat("Diversity and Quantity of Flight Behavior over"," ",num2str(num_days)," ","days of training."));
disp(strcat("Total # Flights:"," ",num2str(num_flights)));
disp(strcat("Total # Unique Transition Types (2-gram):"," ",num2str(number_unique_transitions_2)));
disp(strcat("Total # Unique Transition Types (3-gram):"," ",num2str(number_unique_transitions_3)));
disp(strcat("Total # of 3-gram transition types with >=10 occurances:"," ",num2str(TG_greater_than_10)));
disp(strcat("Total # of 2-gram transition types with >=10 occurances:"," ",num2str(DG_greater_than_10)));
disp(strcat("Total # of 3-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));
disp(strcat("Total # of 2-gram transition types with >=30 occurances:"," ",num2str(TG_greater_than_30)));

clear ROI_Data aligned_data_struct cell_registered_struct master_track_file;
clear TG_greater_than_10 DG_greater_than_10 DG_greater_than_30 TG_greater_than_30;
clear TG_transition_agg DG_transition_agg TG_transition_number_list DG_transition_number_list;
clear CombinedROI ROI_Data FlightPaths c_s_34 FlightAlignedROI;
clear co number_unique_transitions_2 number_unique_transitions_3;

%% Grand totals
t30_sum = grand_totals.ge.subset1.number_unique_TG_GT30 + grand_totals.ga.subset1.number_unique_TG_GT30 + grand_totals.go.subset1.number_unique_TG_GT30 + grand_totals.gi.subset1.number_unique_TG_GT30 + grand_totals.za.subset1.number_unique_TG_GT30 + grand_totals.zu.subset1.number_unique_TG_GT30 + grand_totals.z2.subset1.number_unique_TG_GT30
t10_sum = grand_totals.ge.subset1.number_unique_TG_GT10 + grand_totals.ga.subset1.number_unique_TG_GT10 + grand_totals.go.subset1.number_unique_TG_GT10 + grand_totals.gi.subset1.number_unique_TG_GT10 + grand_totals.za.subset1.number_unique_TG_GT10 + grand_totals.zu.subset1.number_unique_TG_GT10 + grand_totals.z2.subset1.number_unique_TG_GT10
disp(strcat("Total number of unique 3-gram transitions with at least 10 ocurrances:"," ",num2str(t10_sum)));
disp(strcat("Total number of unique 3-gram transitions with at least 30 ocurrances:"," ",num2str(t30_sum)));
disp(strcat("Ge:"," ",num2str(grand_totals.ge.subset1.number_unique_TG_GT30)));
disp(strcat("Ga:"," ",num2str(grand_totals.ga.subset1.number_unique_TG_GT30)));
disp(strcat("Gi:"," ",num2str(grand_totals.gi.subset1.number_unique_TG_GT30)));
disp(strcat("Go:"," ",num2str(grand_totals.go.subset1.number_unique_TG_GT30)));
disp(strcat("Zack:"," ",num2str(grand_totals.za.subset1.number_unique_TG_GT30)));
disp(strcat("Zuzu:"," ",num2str(grand_totals.zu.subset1.number_unique_TG_GT30)));
disp(strcat("Zoom2:"," ",num2str(grand_totals.z2.subset1.number_unique_TG_GT30)));



