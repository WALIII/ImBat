function [ out_data stats] = Imbat_ROI2ManyFlights(FlightAlignedROI);
% Imbat_ROI2ManyFlights




% 01/21/2021
% WAL3


%% To Do:
% - organize the code within 'ImBat_ROI_Behav_Correlation.m'
% -

%% ---------------------------------------------------------------------- %

% User inputs
clusters_2_use = [1 2 3];

% get behavioral correlation for each flight cluster
counter = 1; 
for i = clusters_2_use;
    FlightAlignedROI_combined{1} = FlightAlignedROI{i};
    out{counter} = ImBat_ROI_Behav_Correlation(FlightAlignedROI_combined);
    counter = counter+1;
end

%ImBat_ClusterCalciumVar(FlightAlignedROI,8);
% check if out is empty:
counter = 1;
for i = 1:size(out,2)
    if isempty(out{i})
    else
        out2{counter} = out{i};
        counter = counter+1;
    end
end
try
out = out2;
catch
    out_data = [];
    disp(' no data here...')
    return
end

%% Cocncatonate data:
out_f.D = [];
out_f.calPeaks = [];
out_f.calROI_ID = [];
for iii = 1: size(out,2)
    try
    out_f.D = cat(2,out_f.D,out{iii}.D);
    out_f.calPeaks = cat(2,out_f.calPeaks,out{iii}.calPeaks);
    out_f.calROI_ID = cat(2,out_f.calROI_ID,out{iii}.calROI_ID);
    catch
    end
end

ROI_idx = cell2mat(out_f.calROI_ID);


%% Plot concatonated data
figure();
keep_calPeaks = [];
keep_D = [];
R1 = [];
R2 = [];
hold on;
for i = 1:max(cell2mat(out_f.calROI_ID))
    % first in out1;
    idx2try = find(ROI_idx ==i);
    if size(idx2try,2)>1; % if there is more than one peak for this cell
        counter = 1;
        for ii = 1:size(idx2try,2)
            temp_calPeaks(counter) = mean(((out_f.calPeaks{idx2try(ii)})));
            temp_D(counter) =  median(out_f.D{idx2try(ii)});
            counter = counter+1;
        end
        
         % Make ratio
        temp_R1 = temp_calPeaks(1)./temp_D(1);
        temp_R2 = temp_calPeaks(2)./temp_D(2);
        
        % Subtract One of the Two 
        temp_calPeaks = temp_calPeaks-temp_calPeaks(2);
        temp_D = temp_D-temp_D(2);
       
        plot(temp_D,temp_calPeaks,'ob')
        keep_calPeaks = [keep_calPeaks,temp_calPeaks];
        keep_D = [keep_D,temp_D];
        R1 = [R1, temp_R1];
        R2 = [R2, temp_R2];

        clear temp_D temp_calPeaks a b G1 G2 temp_R1 temp_R2        
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
histogram(G1,'BinWidth',0.25,'Normalization','probability')
histogram(G2,'BinWidth',0.25,'Normalization','probability')
legend('Higher Variance Flights','Lower Varaince flights');
xlabel('     less <---    [ ROI consistancy ]   ---> more ');
ylabel('Frequency');

try
% Stats on the histograms:
[pval_combined_data,~] = ranksum(G1,G2)

out_data.keep_calPeaks = keep_calPeaks;
out_data.keep_D = keep_D;
out_data.G1 = G1;
out_data.G2 = G2;

stats = pval_combined_data;
catch
    
out_data.keep_calPeaks = [];
out_data.keep_D = [];
out_data.G1 = [];
out_data.G2 = [];
stats = [];
end
    