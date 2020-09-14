day = 9;

ROIs_gal = [28 20 1 23 12 22 10 8 11 24 NaN 2 21 30 19;
    3 2 10 28 11 1 5 33 8 35 NaN 6 22 32 29;
    4 5 11 24 5 1 16 10 2 18 14 8 25 19 9;
    11 22 4 18 3 1 14 5 19 39 9 17 36 25 8;
    14 3 16 21 2 1 5 7 8 26 NaN 9 27 6 4;
    5 13 41 23 1 21 3 24 6 22 2 25 16 15 7;
    12 3 34 19 2 14 6 15 9 36 5 10 35 20 1;
    25 26 16 32 1 12 4 19 5 28 15 NaN 34 3 2;
    32 34 29 51 7 10 6 40 16 45 5 8 42 26 43];  % 15 stable manually selected ROIs across 9 days for Gal

for ROI_i = 1:length(ROIs_gal(day,:))
    ROIdir = dir(['*ROI_' num2str(ROIs_gal(day,ROI_i)) '_cluster*']);
    for copy_i = 1:length(ROIdir)
        copyfile(ROIdir(copy_i).name,'X:\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\200909-placeCells across days\');
    end
    
end
