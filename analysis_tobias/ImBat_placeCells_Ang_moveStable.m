nDays = 1;
batId = 'Gen';

if strcmp(batId,'Gal') 
% 15 stable manually selected ROIs across 9 days for Gal
ROIs_manual = [28 20 1 23 12 22 10 8 11 24 NaN 2 21 30 19;
    3 2 10 28 11 1 5 33 8 35 NaN 6 22 32 29;
    4 5 11 24 5 1 16 10 2 18 14 8 25 19 9;
    11 22 4 18 3 1 14 5 19 39 9 17 36 25 8;
    14 3 16 21 2 1 5 7 8 26 NaN 9 27 6 4;
    5 13 41 23 1 21 3 24 6 22 2 25 16 15 7;
    12 3 34 19 2 14 6 15 9 36 5 10 35 20 1;
    25 26 16 32 1 12 4 19 5 28 15 NaN 34 3 2;
    32 34 29 51 7 10 6 40 16 45 5 8 42 26 43]; 
g = dir('Ga*');
elseif strcmp(batId,'Gen') 
% 20 stable manually selected ROIs across 5 days for Gen
ROIs_manual = [NaN NaN 10 3 16 12 17 18 27 29 8 9 NaN NaN 21 11 31 15 20 25;
    8 17 5 1 2 6 21 10 18 31 NaN 11 51 53 28 4 38 19 23 20;
    50 54 12 3 48 18 27 15 31 34 NaN NaN 28 NaN 29 25 24 22 38 14;
    8 NaN 4 28 3 18 10 35 42 25 13 NaN 50 39 46 NaN 49 2 32 26;
    14 NaN 3 28 2 6 33 26 18 45 NaN NaN 25 NaN 32 NaN 37 8 28 11];
g = dir('Ge*');
end
% ROIs_gal = [28 20 1 23 12 22 10 8 11 24 NaN 2 21 30 19;
%     3 2 10 28 11 1 5 33 8 35 NaN 6 22 32 29;
%     4 5 11 24 5 1 16 10 2 18 14 8 25 19 9;
%     11 22 4 18 3 1 14 5 19 39 9 17 36 25 8;
%     14 3 16 21 2 1 5 7 8 26 NaN 9 27 6 4;
%     5 13 41 23 1 21 3 24 6 22 2 25 16 15 7;
%     12 3 34 19 2 14 6 15 9 36 5 10 35 20 1;
%     25 26 16 32 1 12 4 19 5 28 15 NaN 34 3 2;
%     32 34 29 51 7 10 6 40 16 45 5 8 42 26 43];  % 15 stable manually selected ROIs across 9 days for Gal
for day_i = nDays
for ROI_i = 1:length(ROIs_manual(day_i,:))
    ROIdir = dir(['*ROI_' num2str(ROIs_manual(nDays,ROI_i)) '_cluster*']);
    for copy_i = 1:length(ROIdir)
        copyfile(ROIdir(copy_i).name,'X:\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\201001\placeCells across days\');
    end
    
end
end
