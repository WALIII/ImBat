function ImBat_PlotRedBlue(out_data)

% plot the red, blue stability plots seperately and togetehr..

% WAL3
% d11/30/2020


% out_data = ImBat_analysis_20201029(CombinedROI,flightPaths,2)

%[out_final] = ImBat_analysis_11122020(CombinedROI,flightPaths);
% Then, combine data:


figure(); 
hold on;
ROI_ON = out_data.ROI_ON;
% get bound to plot
bound2plot = 1:500;

col = [1,0,0; 0,0,1];
counter = 1; % reset counter
for i = 1:size(out_data.x_2save,1);
    

    for ii = 1
        x1 = out_data.x_2save(i,ii);
         y1 = i;%y1 = out_data.y_2save(i,ii);
        err1 = out_data.err_2save(i,ii);
        errorbar(x1,y1,err1,'horizontal','o','color',col(ii,:),'CapSize',1);
        
    end
    counter = counter+1;
    
end
xlim([ 0 bound2plot(end)]);
ylim([0 i+1]);

% Set where ticks will be
% Get axis handle
ax = gca;
ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
% Set TickLabels;
ylabel('filghts')
xlabel('time from takeoff');
ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};



plot([ROI_ON ROI_ON],[0 i+1],'--k')
%plot([ROI_ON+out_data.fl_length ROI_ON+out_data.fl_length],[0 i+1],'--k')



figure(); 
hold on;
ROI_ON = out_data.ROI_ON;
% get bound to plot
bound2plot = 1:500;

col = [1,0,0; 0,0,1];
counter = 1; % reset counter
for i = 1:size(out_data.x_2save,1);
    

    for ii = 2
        x1 = out_data.x_2save(i,ii);
         y1 = i;%y1 = out_data.y_2save(i,ii);
        err1 = out_data.err_2save(i,ii);
        errorbar(x1,y1,err1,'horizontal','o','color',col(ii,:),'CapSize',1);
        
    end
    counter = counter+1;
    
end
xlim([ 0 bound2plot(end)]);
ylim([0 i+1]);

% Set where ticks will be
% Get axis handle
ax = gca;
ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
% Set TickLabels;
ylabel('filghts')
xlabel('time from takeoff');
ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};



plot([ROI_ON ROI_ON],[0 i+1],'--k')
%plot([ROI_ON+out_data.fl_length ROI_ON+out_data.fl_length],[0 i+1],'--k')




figure();
hold on;

counter = 1; % reset counter
for i = 1:size(out_data.x_2save,1);
    

    for ii = 1:2
        x1 = out_data.x_2save(i,ii);
        y1 = i; %y1 = out_data.y_2save(i,ii);
        err1 = out_data.err_2save(i,ii);
        errorbar(x1,y1,err1,'horizontal','o','color',col(ii,:),'CapSize',1);
        
    end
    counter = counter+1;
    
end
xlim([ 0 bound2plot(end)]);
ylim([0 i+1]);

% Set where ticks will be
% Get axis handle
ax = gca;
ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
% Set TickLabels;
ylabel('filghts')
xlabel('time from takeoff');
ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};



plot([ROI_ON ROI_ON],[0 i+1],'--k')
%plot([ROI_ON+out_data.fl_length ROI_ON+out_data.fl_length],[0 i+1],'--k')


figure(); 
histogram(cat(2,out_data.x_2save(:,1),out_data.x_2save(:,2)),'binwidth',15);

% Stats
X1 = abs((out_data.x_2save(:,1)-out_data.x_2save(:,2))/30);
% g = find(X1<2);
% stable_per = (size(g,1)./size(X1,1))*100;

X2 = out_data.err_2save(:,1)+out_data.err_2save(:,2);

stable_per_error_X = find((X1-X2/2)<0);
stable_per_error = (size(stable_per_error_X,1)./size(X1,1))*100


