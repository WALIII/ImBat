function Imbat_ROI2ManyFlights(FlightAlignedROI);


% for analyis perfomed on Jan 20th, 2021
% WAL3

clusters_2_use = 3

for i = 1:clusters_2_use;
FlightAlignedROI_combined{1} = FlightAlignedROI{i};
out{i} = ImBat_ROI_Behav_Correlation(FlightAlignedROI_combined);
end

%ImBat_ClusterCalciumVar(FlightAlignedROI,8);


% cocncat data:
out_f.D = [];
out_f.calPeaks = [];
out_f.calROI_ID = [];
for iii = 1: size(out,2)
   out_f.D = cat(2,out_f.D,out{iii}.D);
   out_f.calPeaks = cat(2,out_f.calPeaks,out{iii}.calPeaks);
   out_f.calROI_ID = cat(2,out_f.calROI_ID,out{iii}.calROI_ID);
end

ROI_idx = cell2mat(out_f.calROI_ID);

figure(); 
keep_calPeaks = [];
keep_D = [];
hold on;
for i = 1:max(cell2mat(out_f.calROI_ID))
    % first in out1;
    idx2try = find(ROI_idx ==i)
    if size(idx2try,2)>1;
        counter = 1;
        for ii = 1:size(idx2try,2)
         temp_calPeaks(counter) = mean(((out_f.calPeaks{idx2try(ii)})));
         temp_D(counter) =  median(out_f.D{idx2try(ii)});
         counter = counter+1;
        end
        
                % subtract the smallest
%         [a b] = min(temp_calPeaks);
        temp_calPeaks = temp_calPeaks-temp_calPeaks(2);
        temp_D = temp_D-temp_D(2);
        plot(temp_D,temp_calPeaks,'o')
        keep_calPeaks = [keep_calPeaks,temp_calPeaks];
        keep_D = [keep_D,temp_D];
            clear temp_D temp_calPeaks a b G1 G2
        % subtract the smallest 
        
        
    end
end
plot([0 0],[-2 2],'--r')  
plot([-600 600],[0 0],'--r')
xlabel('relative change in mean flight variance from center at peak')
ylabel('relative change in df/f');

figure(); 
hold on;
G1 = keep_calPeaks(find(keep_D>0));
G2 = keep_calPeaks(find(keep_D<0));
[pval_combined_data,~] = ranksum(G1,G2)
histogram(G1,'BinWidth',0.25,'Normalization','probability')
histogram(G2,'BinWidth',0.25,'Normalization','probability')
legend('Higher Variance Flights','Lower Varaince flights');
xlabel('     less <---    [ ROI consistancy ]   ---> more ');
ylabel('Frequency');


    
    