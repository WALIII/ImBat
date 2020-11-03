function ImBat_analysis_10212020(flightPaths,ROI_Data,CombinedROI,clst);


%clst = 2;

[FlightAlignedROI] = ImBat_Align_FC(CombinedROI,flightPaths,clst);


figure();
hold on;

stop_time = 1;
sub2plot = 5;


% Plot all flights for 5 days, hilighting the clustered flights
A = flightPaths.tracjectoriesRaw*1000;

for i = 1:sub2plot%max(flightPaths.day); % number of days
    subplot(1,sub2plot,i);
    start_time = stop_time;
    stop_time = start_time+length(ROI_Data{i}.Alignment.out.flights(:,1))-1;
    % plot all flights for one day:
    bound = start_time:stop_time;
    
    hold on;
    % plot all flights
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'k');
    plot1.Color(4) = 0.4;
    % plot clust flights
    dayFlights = flightPaths.day(flightPaths.clusterIndex{clst});
    todaysFlight = find(dayFlights==i);
    for ii = 1:size(todaysFlight,1)
        plot3(squeeze(FlightAlignedROI.ClustFlight(:,1,todaysFlight(ii))),squeeze(FlightAlignedROI.ClustFlight(:,2,todaysFlight(ii))),squeeze(FlightAlignedROI.ClustFlight(:,3,todaysFlight(ii))),'LineWidth',2,'color','r')
    end
    grid on;
    view( -37.5000,30)
    
end

% Stats on the clustered fights


% %% ROI Analysis
CutCells = FlightAlignedROI.C_raw;
ROI_ON = FlightAlignedROI.ROI_ON;

% get dates:
transition_points = find((diff(FlightAlignedROI.CutCells_date)>=1));  %?
transition_points = [1 transition_points size(CutCells,3)];


% get bound to plot
bound2plot = 1:400;
counter = 1;
col = hsv(size(transition_points,2)+1);
figure();
for i = 1:size(CutCells,1);
    subplot(10,1,1:2);% plot shadding for each cell, over lay for each day
hold on;
    for ii = 1:size(transition_points,2)-1
        
        adata = zscore(squeeze(CutCells(i,bound2plot,transition_points(ii):transition_points(ii+1))),[],1)';
        
        L = size(adata,2);
        se = std(adata)/2;%/10;%sqrt(length(adata));
        mn = nanmean(adata);
        h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col(ii,:)); alpha(0.5);
        plot(mn,'Color',col(ii,:));
    end
    hold off
    subplot(10,1,3:10);
    
    adata = zscore(squeeze(CutCells(i,bound2plot,:)),[],1)';
    imagesc(adata);
    % draw date lines:
    hold on;
    for ii = 1:size(transition_points,2)
        line([1,size(CutCells,2)], [transition_points(ii),transition_points(ii)], 'Color','r','LineWidth',2);
    end
    % Get axis handle
    ax = gca;
    % Set where ticks will be
    ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
    % Set TickLabels;
    ylabel('filghts')
    xlabel('time from takeoff');
    ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};
    title([ 'ROI ' num2str(i)]);
    pause();
    clf('reset');
    clear adata
end





% %% ROI Analysis
% col = hsv(6);
% CutCells = FlightAlignedROI.C;
% ROI_ON = FlightAlignedROI.ROI_ON;
% figure();
% hold on;
% counter = 1;
% for i = [1 3 8 10 11 15 17 20 22 32 37 38]
%     % draw date lines:
%     for ii = 1:size(transition_points,2)-1
%         
%         adata = zscore(squeeze(CutCells(i,bound2plot,transition_points(ii):transition_points(ii+1))),[],1)'+counter*4;
%         
%         L = size(adata,2);
%         se = std(adata)/2;%/10;%sqrt(length(adata));
%         mn = nanmean(adata);
%         h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col(ii,:)); alpha(0.5);
%         plot(mn,'Color',col(ii,:));
%     end
%     counter = counter+1;
% end




