
% Light Dark Analysis 

% Load in Data from Light/Dark Session:

L = [];
D = [];
R = [];
for i = 1:3
    try
    out_Ge{i} = ImBat_PlotAlignedROIs(FlightAlignedROI{i},ROI_Data,flightPaths,'sort',dark_cluster{i+1},'smooth',1,'display',0); % cluster 2 is # 2 in FLightAlighedROI

    [Ltemp, Dtemp, Rtemp] = Scrap_remap(out_Ge{i});
L = cat(2,L,Ltemp);
D = cat(2,D,Dtemp);
R = cat(2,R,Rtemp);
    catch
        disp('no flights');
    end
    
end


figure(); 
hold on;
histogram(1-R,'BinWidth',0.10,'Normalization','probability','FaceColor','r');
histogram(1-D,'BinWidth',0.10,'Normalization','probability','FaceColor','k');
histogram(1-L,'BinWidth',0.10,'Normalization','probability','FaceColor','c');

[pval_combined_data,~] = ranksum(L,D);
[pval_combined_data_LR,~] = ranksum(L,R);

[pval_combined_data_DR,~] = ranksum(D,R);

disp(['pval of Light v Dark = ', num2str(pval_combined_data)]);
disp(['pval of Light v Shiffle = ', num2str(pval_combined_data_LR)]);
disp(['pval of Dark v Shiffle = ', num2str(pval_combined_data_DR)]);


legend('Shuffled ID','Dark vs Light', 'Light vs Light''')
title( 'Remapping Index [ 1-r1_d ]')
xlabel('Remapping index');
ylabel('% ROIs');
