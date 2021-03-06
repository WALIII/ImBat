function [out_final] = ImBat_analysis_11122020(CombinedROI,flightPaths);
% Generate source data for Red-blue plots

% WAL3
% 3/16/2021

clear out_1 out_2 out_3 out3
out_1 = ImBat_analysis_20201029(CombinedROI,flightPaths,2);
out_2 = ImBat_analysis_20201029(CombinedROI,flightPaths,3);
out_3 = ImBat_analysis_20201029(CombinedROI,flightPaths,4);
close all;
% combine data


        out3.x1 = cat(1,out_1.x_2save,out_2.x_2save,out_3.x_2save);
        out3.y1 = cat(1,out_1.y_2save,out_2.y_2save,out_3.y_2save);
        out3.err1 = cat(1,out_1.err_2save,out_2.err_2save,out_3.err_2save);
        out3.FL = cat(1,ones(size(out_1.err_2save,1),1)*out_1.fl_length, ones(size(out_2.err_2save,1),1)*out_2.fl_length,ones(size(out_3.err_2save,1),1)*out_3.fl_length);
        out_final.x_2save = out3.x1;
        out_final.y_2save = out3.y1;
        out_final.err_2save = out3.err1;
        out_final.FL = out3.FL;
        out_final.ROI_ON = out_1.ROI_ON;
        
% resort based on max location:
        

[~, sort_idx] = sort(out3.x1(:,1));
out3.x1 = out3.x1(sort_idx,:);
out3.y1 = out3.y1 (sort_idx,:);
out3.err1 = out3.err1(sort_idx,:);
out3.FL = out3.FL(sort_idx,:);

% EXPORT:
        out_final.forward.x_2save = out3.x1;
        out_final.forward.y_2save = out3.y1;
        out_final.forward.err_2save = out3.err1;
        out_final.forward.FL = out3.FL;
        out_final.forward.ROI_ON = out_1.ROI_ON;
        
% Plot data
col = [1,0,0; 0,0,1];
ROI_ON = out_1.ROI_ON;

figure(); 
histogram(cat(2,out3.x1(:,1),out3.x1(:,2)),'binwidth',15);
hold on;
for i = 1:size(out3.x1,1) % for all pairs:
        for ii = 1:2
        x1 = out3.x1(i,ii);
        y1 = i;
        err1 = out3.err1(i,ii);
        errorbar(x1,y1,err1,'horizontal','.','color',col(ii,:),'CapSize',1,'MarkerSize',10);

        end
end
ax = gca;
ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
% Set TickLabels;
ylabel('ROIs')
xlabel('time from takeoff');
ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};





%% now plot by the reward

out3.x1 = out3.x1-out3.FL+1;

[~, sort_idx] = sort(out3.x1(:,1),'descend');
out3.x1 = out3.x1(sort_idx,:);
out3.y1 = out3.y1 (sort_idx,:);
out3.err1 = out3.err1(sort_idx,:);
out3.FL = out3.FL(sort_idx,:);

% EXPORT:
        out_final.reverse.x_2save = out3.x1;
        out_final.reverse.y_2save = out3.y1;
        out_final.reverse.err_2save = out3.err1;
        out_final.reverse.FL = out3.FL;
        out_final.reverse.ROI_ON = out_1.ROI_ON;
% Plot data
col = [1,0,0; 0,0,1];
ROI_ON = out_1.ROI_ON;

figure(); 
histogram(cat(2,out3.x1(:,1),out3.x1(:,2)),'binwidth',15);
hold on;
for i = 1:size(out3.x1,1) % for all pairs:
        for ii = 1:2
        x1 = out3.x1(i,ii);
        y1 = i;
        err1 = out3.err1(i,ii);
        errorbar(x1,y1,err1,'horizontal','.','color',col(ii,:),'CapSize',1,'MarkerSize',10);

        end
end
ax = gca;
ax.XTick = [ROI_ON-150 ROI_ON-120 ROI_ON-90 ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60 ROI_ON+90 ROI_ON+120];
% Set TickLabels;
ylabel('ROIs')
xlabel('time from Landing');
ax.XTickLabel = {'-5','-4','-3','-2','-1','0','1','2','3','4'};




%% 
% plot difference by initial field location:
figure(); 
plot((out3.x1(:,1)-out_1.ROI_ON)/30,abs(out3.x1(:,1)-out3.x1(:,2))/30,'*');
ylabel('Seconds of change');
xlabel('time of peak relative to flight');


figure(); 
plot((out3.err1(:,1)+out3.err1(:,2))/30,abs(out3.x1(:,1)-out3.x1(:,2))/30,'*');
xlabel('Summed Error ');
ylabel('Seconds of change');


% EXPORT DATA:





