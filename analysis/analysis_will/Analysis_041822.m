
% Look at the 'stabilization' of place cells over days

% get 'structured' data:
[FlightAlignedROI_Clustered] = ImBat_Combine_Clusters(CombinedROI,flightPaths,2:10);

% get 'unstructured' data:
FlightAlignedROI_Unclustered{1} = ImBat_Align_FC(CombinedROI,flightPaths,1);

% combine bothL
FlightAlignedROI_All{1}.S = cat(3,FlightAlignedROI_Clustered{1}.S, FlightAlignedROI_Unclustered{1}.S);
FlightAlignedROI_All{1}.ClustFlight_withPads  = cat(3,FlightAlignedROI_Clustered{1}.ClustFlight_withPads,FlightAlignedROI_Unclustered{1}.ClustFlight_withPads); 
FlightAlignedROI_All{1}.cluster_idX  = cat(2,FlightAlignedROI_Clustered{1}.cluster_idX,FlightAlignedROI_Unclustered{1}.cluster_idX); 
FlightAlignedROI_All{1}.cluster_ID  =  cat(2,FlightAlignedROI_Clustered{1}.cluster_ID,ones(1,size(FlightAlignedROI_Unclustered{1}.ClustFlight_withPads,3))*FlightAlignedROI_Unclustered{1}.clust_number);
FlightAlignedROI_All{1}.CutCells_date  =  cat(2,FlightAlignedROI_Clustered{1}.CutCells_date,FlightAlignedROI_Unclustered{1}.CutCells_date);


% plot example cells
warning off;
clear val2use val2use_rand
for ROI2use = 1:100
[R,Map2save,occupancyMap2save,days2use] = ImBat_ratemeaps(FlightAlignedROI_Clustered,'roi2plot',ROI2use,'cluster',1,'randomize',0);
val2use{ROI2use} = squeeze(max(max(Map2save)));
clear Map2save
% randomize
for itter = 1:10;
    itter
[R,Map2save,occupancyMap2save,days2use] = ImBat_ratemeaps(FlightAlignedROI_Clustered,'roi2plot',ROI2use,'cluster',1,'randomize',1);
val2use_rand{ROI2use}(:,itter) = squeeze(max(max(Map2save)));
clear Map2save
end
figure();
hold on;
plot(val2use{ROI2use},'b');
plot(val2use_rand{ROI2use},'r');
end

% calculate significance

for ROI2use = 1:10
    g = val2use{ROI2use};
    g2 = val2use_rand{ROI2use}';
   
   for day2use = 1:35
       if isnan(g(day2use))
        vals(day2use,ROI2use) = NaN; 
       else
  vals(day2use,ROI2use) =  size(find(g2(day2use,:)>g(day2use)),2)./size(g2(day2use,:),2);
       end
   end
end

% plot percentages
for day2use = 1:35;
prtg(day2use) = size(find(vals(day2use,:)<0.05),2)./(size(find(vals(day2use,:)<0.05),2)+size(find(vals(day2use,:)>0.05),2))
end

figure(); bar(prtg);
title('percentage of sig tuned cells');


   