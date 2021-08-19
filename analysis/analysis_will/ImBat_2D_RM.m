function ImBat_2D_RM;
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
   
disp('Loading data...');
load([subFolders(ii).name, '/CellReg_files/ROI_Data/Saved_Data/Aligned_Data.mat'],'CombinedROI','flightPaths','FlightAlignedROI')



[FlightAlignedROI_rando] = ImBat_Align_FC(CombinedROI,flightPaths,1);
FlightAlignedROI_rando1{1} = FlightAlignedROI_rando;
for i = 1:size(FlightAlignedROI{1}.C,1)-10;
[R_1{ii}(i,:)] = ImBat_Compare2D_ratemeaps_acrossDays(FlightAlignedROI,i,1);
[R_2{ii}(i,:)] = ImBat_Compare2D_ratemeaps_acrossDays(FlightAlignedROI_rando1,i,1);
close all
end

clear CombinedROI flightPaths FlightAlignedROI FlightAlignedROI_rando1
end

% concat together:
Rr = zeros(1,20); Rr(Rr ==0) = NaN; % make empty nan box
Rr2 = Rr;


 figure();
 hold on;
for i = 1:8
try
R = R_1{i};
R2 = R_2{i};
R2(R <.3) = NaN;
R2(R > .95) = NaN;
R(R <.3) = NaN;
R(R > .95) = NaN;

idx = find(R(:,1)>0.01 );
y1 = nanmean(R(idx,:));
x1 = 1:size(R(idx,:),2);
err1 = nanstd(R(idx,:))/2;
errorbar(x1,y1,err1,'LineWidth',2,'color','r');
hold on;

R2(:,1) = (R2(:,1)+ R(:,1))/2;
y2 = nanmean(R2(idx,:));
x2 = 1:size(R2(idx,:),2);
err2 = nanstd(R2(idx,:))/2;
errorbar(x2,y2,err2,'LineWidth',2,'color','b');
clear R R2 x1 x2 y1 y2 err1 err2
catch
   disp('missing data for this bat'); 
end
end
title('Change in 2D corr given stereotyped (r), or unique (b) flights')

ylabel('corr to day 1 (r)');
xlabel('days');

