function [out, stats,meanF,meanC,stats2keep] = ImBat_ROI2Behav(FlightAlignedROI,cell_registered_struct);
% Breakout for computing the ROI to behavior correlation for our dataset

stats = [];
out = [];
meanF = [];
meanC = [];
stats2keep = [];
remove_trans = 0;


% I.  first, get the real ROIs for this dataset, and  prune the
% flightAligned ROIs of ROIs that only exist on one day..
for iii = 1:3;
if remove_trans ==1;
FlightAlignedROI2 = FlightAlignedROI;
G = cell_registered_struct.cell_to_index_map;
G(G>0) = 1;
G2 = sum(G');
% G2(G2>1) = [];
for i = 1:size(FlightAlignedROI,2)
FlightAlignedROI{i}.C_raw = FlightAlignedROI{i}.C_raw(:,:,G2(G2>1));
FlightAlignedROI{i}.C = FlightAlignedROI{i}.C(:,:,G2(G2>1));
FlightAlignedROI{i}.S = FlightAlignedROI{i}.S(:,:,G2(G2>1));
FlightAlignedROI{i}.ClustFlight_withPads = FlightAlignedROI{i}.ClustFlight_withPads(:,:,G2(G2>1));
end
end

dat2use{1} = FlightAlignedROI{iii};
[out2] = ImBat_ROI_Behav_Correlation(dat2use);

meanF = cat(2,meanF,out2.meanF);
meanC = cat(2,meanC,out2.meanC);
stats2keep = cat(2,stats2keep,out2.stats2keep);

out{iii} = out2;
end


% stats
% ========== stats ( regression) ===========%
% clear y x1 X
% X1 = meanF';
% y = mat2gray(meanC');
% x1 = ones(size(X1,1),1);
% X = [x1 X1];    % Includes column of ones
% 
% [~,~,~,~,stats] = regress(y,X);
% stats
% X = [X1];
% mdl = fitlm(X,y);
% figure();
% plotAdded(mdl)
% title('Calcium consistancy vs Flight Variance');
% ylabel(' Mean, Normalized Ca+ activity');
% xlabel('mean Flight distance (cm);');







% II. OPTIONAL: find cells with high probability tracking transitions 


% III. run the correlation for each dataset for the top flightpaths 



% IV. aggregate the stats


