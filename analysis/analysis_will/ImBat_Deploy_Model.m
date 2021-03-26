function [output] = ImBat_Deploy_Model
% Temp function to aggregate and quantify ROI stability over days

% WAL3
% 02/28/2021

% nb: Run in data in dir with 'processed' folders...


%% Aggregate Data

DIR = pwd;
files = dir(pwd);
files(ismember( {files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.

for ii = 1 : length(subFolders)
   
disp('Loading data');
load([subFolders(ii).name, '/CellReg_files/ROI_Data/Saved_Data/Aligned_Data.mat'],'CombinedROI','flightPaths')

for i = 1:max(flightPaths.day);
    
[R, R3, M] = ImBat_Flight_Model_AcrossDays(flightPaths,i);
% mR(i) = nanmean(R);
% stdR(i) = nanstd(R);
% mR3(i) = nanmean(R3);
% stdR3(i) = nanstd(R3);
% store data:
R_store{ii}{i} = R;
R3_store{ii}{i} = R3;
clear R R3 M
end

% shuffle day index
flightPaths2 = flightPaths;
dc1 = flightPaths.day;

% % randomize index
fakeInd = randperm((size(dc1,1)));
dc1 = dc1(fakeInd);
flightPaths2.day = dc1;
for i = 1:max(flightPaths.day);
    [R_s, R3_s, M] = ImBat_Flight_Model_AcrossDays(flightPaths2,i);
% mR2(i) = nanmean(R);
% stdR2(i) = nanstd(R);
% mR32(i) = nanmean(R3);
% stdR32(i) = nanstd(R3);

% store data:
R_store_s{ii}{i} = R_s;
R3_store_s{ii}{i} = R3_s;
clear R2 R32 M2
end
end

% Export data
output.R_store = R_store;
output.R3_store = R3_store;
output.R_store_s = R_store_s;

% combine data
clear McR StdCR 
for ii = 1:12 
    c_R = []; c_RS = []; c_RShuff = [];

    for i = 1:size(R_store_s,2)
   try
   c_R = cat(2,c_R,R_store{i}{ii}); % simulated data
   c_RS = cat(2,c_RS,R_store_s{i}{ii}); % sorted data
   c_RShuff = cat(2,c_RShuff,R3_store_s{i}{ii}); % shuffeled data

    catch
    c_R = cat(2,c_R,NaN);
    end
    end
    % combine simulated data
McR(ii) = nanmean(c_R);
StdCR(ii) = nanstd(c_R);
% combine flight shuffle data
McR_S(ii) = nanmean(c_RS);
StdCR_S(ii) = nanstd(c_RS);
% combine ROI shuffle data
McR_Shuff(ii) = nanmean(c_RShuff);
StdCR_Shuff(ii) = nanstd(c_RShuff);
clear c_R
    end 
figure();
hold on;
errorbar(1:length(McR),McR,StdCR/4,'LineWidth',2);
errorbar(1:length(McR_Shuff),McR_Shuff,StdCR_Shuff/4,'LineWidth',2);
errorbar(1:length(McR_S),McR_S,StdCR_S/4,'LineWidth',2);


ylim([-0.4 1.1]);
xlim([0.9 12.2]);
xlabel('days');
ylabel('corr to day 1');
title('Apparent Stability Simulation'); 
legend('simulated data','shuffled ROIs','shuffled Flights');