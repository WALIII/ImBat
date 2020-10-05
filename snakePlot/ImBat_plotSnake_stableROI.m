%function to plot snakePlots flight paths of the
%bats for each day focusing only on the stable neurons from ROIs_manual
batId = 'Gen';
clustNum = 3;
saveFlag = 1; %do you want to save the figures and output structure?
saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
% Check if folder exists
if exist([saveDir1 datestr(now,'yymmdd') filesep 'snakePlots'])>0;
    disp('Youve been working today..');
else
    mkdir([saveDir1 datestr(now,'yymmdd') filesep 'snakePlots'])
end
saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'snakePlots' '\'];

xlimFlight = 1300;
xlimCalcium = xlimFlight/4;
xlimFlightPre = 800;
xlimCalciumPre = xlimFlightPre/4;
xlimFlightPost = 800;
xlimCalciumPost = xlimFlightPost/4;
if strcmp(batId,'Gal') 
% 15 stable manually selected ROIs across 9 days for Gal
ROIs_manual = [28 20 1 23 12 22 10 8 11 24 NaN 2 21 30 19;
    3 2 10 28 11 1 5 33 8 35 NaN 6 22 32 29;
    4 5 11 24 5 1 16 10 2 18 14 8 25 19 9;
    11 22 4 18 3 1 14 5 19 39 9 17 36 25 8;
    14 3 16 21 2 1 5 7 8 26 NaN 9 27 6 4;
    5 13 41 23 1 21 3 24 6 22 2 25 16 15 7;
    12 3 34 19 2 14 6 15 9 36 5 10 35 20 1;
    25 26 16 32 1 12 4 19 5 28 15 NaN 34 3 2;
    32 34 29 51 7 10 6 40 16 45 5 8 42 26 43]; 
g = dir('Ga*');
elseif strcmp(batId,'Gen') 
% 20 stable manually selected ROIs across 5 days for Gen
ROIs_manual = [NaN NaN 10 3 16 12 17 18 27 29 8 9 NaN NaN 21 11 31 15 20 25;
    8 17 5 1 2 6 21 10 18 31 NaN 11 51 53 28 4 38 19 2 20;
    50 54 12 3 48 18 27 15 31 34 NaN NaN 28 NaN 29 25 24 22 38 14;
    8 NaN 4 28 3 18 10 35 42 25 13 NaN 50 39 46 NaN 49 2 32 26;
    14 NaN 3 28 2 6 33 26 18 45 NaN NaN 25 NaN 32 NaN 37 8 28 11];
g = dir('Ge*');
end
z = dir('Z1*');
dirTop = vertcat(g,z); %find all folders in top quality directory

for d = 1:length(dirTop)
    %load results data
    %extract metadata names and enter processed folder
    cd([dirTop(d).name filesep 'extracted'])
    flyFolders = dir('*fly*extraction');
    batName{d} = flyFolders(end).name(1:3);
    dateSesh{d} = flyFolders(end).name(5:10);
    sessionType{d} = flyFolders(end).name(12:16);
    
    cd(flyFolders(end).name);
    
    %load snakeTrace data
    %extract metadata names and enter analysis folder
    dirAnalysis = dir('analysis_*');
    cd(dirAnalysis(end).name);
    sp = dir('*snakePlotData.mat');
    snakeTraceT(d) = load(sp(end).name);
    fp = dir('*flightPaths.mat');
    flightPathsF(d) = load(fp(end).name);
    close all;
    cd(dirTop(d).folder)
end

%% plot all stable ROIs across all 9 days sorted by order of ROIs 
plotSnake_acrossDays = figure('units','normalized','outerposition',[0 0 0.5 1]);
sgtitle(['Mean Activity Across 5 Days for Same Flight Path: ' batName{d} ' 20319 to 200324']);
col = jet(length(ROIs_manual(:,1)));
for day_i = 1:length(ROIs_manual(:,1))
    %plot flightpath for cluster #2
    for flight_i = 1:length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})
        subplot(3,length(ROIs_manual(:,1)),day_i);
        plot3(flightPathsF(day_i).flightPaths.pos(1,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(2,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(3,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),...
            '-','LineWidth',1,'Color', col(day_i,:));
        hold on;
        view(0,90);
        xlim([-3 3])
        ylim([-3 3])
        title(['Day ' num2str(day_i)]);
        if day_i == 1
            xlabel('m'); ylabel('m');
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
       
        %plot velocity for cluster #2
        subplot(3,length(ROIs_manual(:,1)), length(ROIs_manual(:,1)) + day_i);
        plot(1:length(snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawFlight{clustNum}(flight_i,1:xlimFlight)),snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawFlight{clustNum}(flight_i,1:xlimFlight));
        hold on;
        if day_i == 1
            ylabel('Velocity (m/s)');
            xlabel('Time (s)');
            yt = get(gca,'YTick');
            xt = get(gca,'XTick');
            set(gca,'xticklabel',[]);
            %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        ylim([0 4.5]);
                title([num2str(length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})) ' flights']);
    end
    
    %plot all ROIs normalized to itself and plotted in 1-n numerical order
    %across all days
    subplot(3,length(ROIs_manual(:,1)), 2*length(ROIs_manual(:,1)) + day_i);
    imagesc(snakeTraceT(day_i).snakeTrace.cRaw.normMeanTraceEachFlight{clustNum}(:,1:xlimCalcium),[1 4.5]);
    colormap(hot);
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    else
        set(gca,'xticklabel',[]);
    end 
    hold off
end

% plot all stable ROIs across all 9 days sorted by their own preference
plotSnake_acrossDays_prefEach = figure('units','normalized','outerposition',[0 0 0.5 1]);
sgtitle('Mean Activity Across 9 Days for Same Flight Path (Pref Each Day): Gal 203011 to 200320');
col = jet(length(ROIs_manual(:,1)));
for day_i = 1:length(ROIs_manual(:,1))
    %plot flightpath for cluster #2
    for flight_i = 1:length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})
        subplot(3,length(ROIs_manual(:,1)),day_i);
        plot3(flightPathsF(day_i).flightPaths.pos(1,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(2,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(3,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),...
            '-','LineWidth',1,'Color', col(day_i,:));
        hold on;
        view(0,90);
        xlim([-3 3])
        ylim([-3 3])
        title(['Day ' num2str(day_i)]);
        if day_i == 1
            xlabel('m'); ylabel('m');
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        
        %plot velocity for cluster #2
        subplot(3,length(ROIs_manual(:,1)), length(ROIs_manual(:,1)) + day_i);
        plot(1:length(snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawFlight{clustNum}(flight_i,1:xlimFlight)),snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawFlight{clustNum}(flight_i,1:xlimFlight));
        hold on;
        if day_i == 1
            ylabel('Velocity (m/s)');
            xlabel('Time (s)');
            yt = get(gca,'YTick');
            xt = get(gca,'XTick');
            set(gca,'xticklabel',[]);
            %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        ylim([0 4.5]);
                title([num2str(length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})) ' flights']);
    end
    
    %plot all ROIs according to preferred sorting individually
    subplot(3,length(ROIs_manual(:,1)), 2*length(ROIs_manual(:,1)) + day_i);
    imagesc(snakeTraceT(day_i).snakeTrace.cRaw.normTraceFlight{clustNum}(snakeTraceT(day_i).snakeTrace.cRaw.IFlight{clustNum},1:xlimCalcium),[1 4.5]);
    colormap(hot);
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    else
        set(gca,'xticklabel',[]);
    end 
    hold off
end

% plot each ROI across all 9 days in same image with different images for each ROI sorted by day 1 preference
plotSnake_acrossDays_pref1 = figure('units','normalized','outerposition',[0 0 0.5 1]);
sgtitle('Mean Activity Across 9 Days for Same Flight Path (Pref Day 1): Gal 203011 to 200320');
col = jet(length(ROIs_manual(:,1)));
for day_i = 1:length(ROIs_manual(:,1))
    %plot flightpath for cluster #2
    for flight_i = 1:length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})
        subplot(3,length(ROIs_manual(:,1)),day_i);
        plot3(flightPathsF(day_i).flightPaths.pos(1,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(2,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(3,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),...
            '-','LineWidth',1,'Color', col(day_i,:));
        hold on;
        view(0,90);
        xlim([-3 3])
        ylim([-3 3])
        title(['Day ' num2str(day_i)]);
        if day_i == 1
            xlabel('m'); ylabel('m');
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        
        %plot velocity for cluster #2
        subplot(3,length(ROIs_manual(:,1)), length(ROIs_manual(:,1)) + day_i);
        plot(1:length(snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawFlight{clustNum}(flight_i,1:xlimFlight)),snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawFlight{clustNum}(flight_i,1:xlimFlight));
        hold on;
        if day_i == 1
            ylabel('Velocity (m/s)');
            xlabel('Time (s)');
            yt = get(gca,'YTick');
            xt = get(gca,'XTick');
            set(gca,'xticklabel',[]);
            %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        ylim([0 4.5]);
                title([num2str(length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})) ' flights']);
    end
    
    %plot all ROIs according to preferred sorting individually
    subplot(3,length(ROIs_manual(:,1)), 2*length(ROIs_manual(:,1)) + day_i);
    imagesc(snakeTraceT(day_i).snakeTrace.cRaw.normTraceFlight{clustNum}(snakeTraceT(1).snakeTrace.cRaw.IFlight{clustNum},1:xlimCalcium),[1 4.5]);
    colormap(hot);
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    else
        set(gca,'xticklabel',[]);
    end 
    hold off
end

%% Pre-flight activity
%plot all stable ROIs across all 9 days sorted by order of ROIs 
plotSnake_acrossDays_pre = figure('units','normalized','outerposition',[0 0 0.5 1]);
sgtitle('Mean Pre-Activity Across 9 Days for Same Flight Path: Gal 203011 to 200320');
col = jet(length(ROIs_manual(:,1)));
for day_i = 1:length(ROIs_manual(:,1))
    %plot flightpath for cluster #2
    for flight_i = 1:length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})
        subplot(3,length(ROIs_manual(:,1)),day_i);
        plot3(flightPathsF(day_i).flightPaths.pos(1,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(2,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(3,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),...
            '-','LineWidth',1,'Color', col(day_i,:));
        hold on;
        view(0,90);
        xlim([-3 3])
        ylim([-3 3])
        title(['Day ' num2str(day_i)]);
        if day_i == 1
            xlabel('m'); ylabel('m');
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
       
        %plot velocity for cluster #2
        subplot(3,length(ROIs_manual(:,1)), length(ROIs_manual(:,1)) + day_i);
        plot(1:length(snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawPre{clustNum}(flight_i,1:xlimFlightPre)),snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawPre{clustNum}(flight_i,1:xlimFlightPre));
        hold on;
        if day_i == 1
            ylabel('Velocity (m/s)');
            xlabel('Time (s)');
            yt = get(gca,'YTick');
            xt = get(gca,'XTick');
            set(gca,'xticklabel',[]);
            %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        ylim([0 4.5]);
        title([num2str(length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})) ' flights']);
    end
    
    %plot all ROIs normalized to itself and plotted in 1-n numerical order
    %across all days
    subplot(3,length(ROIs_manual(:,1)), 2*length(ROIs_manual(:,1)) + day_i);
    imagesc(snakeTraceT(day_i).snakeTrace.cRaw.normMeanTraceEachPre{clustNum}(:,1:xlimCalciumPre),[1 4.5]);
    colormap(hot);
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    else
        set(gca,'xticklabel',[]);
    end 
    hold off
end

% plot all stable ROIs across all 9 days sorted by their own preference
plotSnake_acrossDays_prefEach_pre = figure('units','normalized','outerposition',[0 0 0.5 1]);
sgtitle('Mean Pre-Activity Across 9 Days for Same Flight Path (Pref Each Day): Gal 203011 to 200320');
col = jet(length(ROIs_manual(:,1)));
for day_i = 1:length(ROIs_manual(:,1))
    %plot flightpath for cluster #2
    for flight_i = 1:length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})
        subplot(3,length(ROIs_manual(:,1)),day_i);
        plot3(flightPathsF(day_i).flightPaths.pos(1,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(2,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(3,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),...
            '-','LineWidth',1,'Color', col(day_i,:));
        hold on;
        view(0,90);
        xlim([-3 3])
        ylim([-3 3])
        title(['Day ' num2str(day_i)]);
        if day_i == 1
            xlabel('m'); ylabel('m');
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        
        %plot velocity for cluster #2
        subplot(3,length(ROIs_manual(:,1)), length(ROIs_manual(:,1)) + day_i);
        plot(1:length(snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawPre{clustNum}(flight_i,1:xlimFlightPre)),snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawPre{clustNum}(flight_i,1:xlimFlightPre));
        hold on;
        if day_i == 1
            ylabel('Velocity (m/s)');
            xlabel('Time (s)');
            yt = get(gca,'YTick');
            xt = get(gca,'XTick');
            set(gca,'xticklabel',[]);
            %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        ylim([0 4.5]);
                title([num2str(length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})) ' flights']);
    end
    
    %plot all ROIs according to preferred sorting individually
    subplot(3,length(ROIs_manual(:,1)), 2*length(ROIs_manual(:,1)) + day_i);
    imagesc(snakeTraceT(day_i).snakeTrace.cRaw.normTracePre{clustNum}(snakeTraceT(day_i).snakeTrace.cRaw.IPre{clustNum},1:xlimCalciumPre),[1 4.5]);
    colormap(hot);
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    else
        set(gca,'xticklabel',[]);
    end 
    hold off
end

% plot each ROI across all 9 days in same image with different images for each ROI sorted by day 1 preference
plotSnake_acrossDays_pref1_pre = figure('units','normalized','outerposition',[0 0 0.5 1]);
sgtitle('Mean Pre-Activity Across 9 Days for Same Flight Path (Pref Day 1): Gal 203011 to 200320');
col = jet(length(ROIs_manual(:,1)));
for day_i = 1:length(ROIs_manual(:,1))
    %plot flightpath for cluster #2
    for flight_i = 1:length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})
        subplot(3,length(ROIs_manual(:,1)),day_i);
        plot3(flightPathsF(day_i).flightPaths.pos(1,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(2,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(3,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),...
            '-','LineWidth',1,'Color', col(day_i,:));
        hold on;
        view(0,90);
        xlim([-3 3])
        ylim([-3 3])
        title(['Day ' num2str(day_i)]);
        if day_i == 1
            xlabel('m'); ylabel('m');
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        
        %plot velocity for cluster #2
        subplot(3,length(ROIs_manual(:,1)), length(ROIs_manual(:,1)) + day_i);
        plot(1:length(snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawPre{clustNum}(flight_i,1:xlimFlightPre)),snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawPre{clustNum}(flight_i,1:xlimFlightPre));
        hold on;
        if day_i == 1
            ylabel('Velocity (m/s)');
            xlabel('Time (s)');
            yt = get(gca,'YTick');
            xt = get(gca,'XTick');
            set(gca,'xticklabel',[]);
            %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        ylim([0 4.5]);
                title([num2str(length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})) ' flights']);
    end
    
    %plot all ROIs according to preferred sorting individually
    subplot(3,length(ROIs_manual(:,1)), 2*length(ROIs_manual(:,1)) + day_i);
    imagesc(snakeTraceT(day_i).snakeTrace.cRaw.normTracePre{clustNum}(snakeTraceT(1).snakeTrace.cRaw.IPre{clustNum},1:xlimCalciumPre),[1 4.5]);
    colormap(hot);
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    else
        set(gca,'xticklabel',[]);
    end 
    hold off
end

%% Post-flight activity
%plot all stable ROIs across all 9 days sorted by order of ROIs 
plotSnake_acrossDays_post = figure('units','normalized','outerposition',[0 0 0.5 1]);
sgtitle('Mean Post-Activity Across 9 Days for Same Flight Path: Gal 203011 to 200320');
col = jet(length(ROIs_manual(:,1)));
for day_i = 1:length(ROIs_manual(:,1))
    %plot flightpath for cluster #2
    for flight_i = 1:length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})
        subplot(3,length(ROIs_manual(:,1)),day_i);
        plot3(flightPathsF(day_i).flightPaths.pos(1,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(2,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(3,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),...
            '-','LineWidth',1,'Color', col(day_i,:));
        hold on;
        view(0,90);
        xlim([-3 3])
        ylim([-3 3])
        title(['Day ' num2str(day_i)]);
        if day_i == 1
            xlabel('m'); ylabel('m');
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
       
        %plot velocity for cluster #2
        subplot(3,length(ROIs_manual(:,1)), length(ROIs_manual(:,1)) + day_i);
        plot(1:length(snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawPost{clustNum}(flight_i,1:xlimFlightPost)),snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawPost{clustNum}(flight_i,1:xlimFlightPost));
        hold on;
        if day_i == 1
            ylabel('Velocity (m/s)');
            xlabel('Time (s)');
            yt = get(gca,'YTick');
            xt = get(gca,'XTick');
            set(gca,'xticklabel',[]);
            %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        ylim([0 4.5]);
                title([num2str(length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})) ' flights']);
    end
    
    %plot all ROIs normalized to itself and plotted in 1-n numerical order
    %across all days
    subplot(3,length(ROIs_manual(:,1)), 2*length(ROIs_manual(:,1)) + day_i);
    imagesc(snakeTraceT(day_i).snakeTrace.cRaw.normMeanTraceEachPost{clustNum}(:,:),[1 4.5]); %1:xlimCalciumPost),[1 4.5]);
    colormap(hot);
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    else
        set(gca,'xticklabel',[]);
    end 
    hold off
end

% plot all stable ROIs across all 9 days sorted by their own preference
plotSnake_acrossDays_prefEach_post = figure('units','normalized','outerposition',[0 0 0.5 1]);
sgtitle('Mean Post-Activity Across 9 Days for Same Flight Path (Pref Each Day): Gal 203011 to 200320');
col = jet(length(ROIs_manual(:,1)));
for day_i = 1:length(ROIs_manual(:,1))
    %plot flightpath for cluster #2
    for flight_i = 1:length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})
        subplot(3,length(ROIs_manual(:,1)),day_i);
        plot3(flightPathsF(day_i).flightPaths.pos(1,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(2,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(3,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),...
            '-','LineWidth',1,'Color', col(day_i,:));
        hold on;
        view(0,90);
        xlim([-3 3])
        ylim([-3 3])
        title(['Day ' num2str(day_i)]);
        if day_i == 1
            xlabel('m'); ylabel('m');
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        
        %plot velocity for cluster #2
        subplot(3,length(ROIs_manual(:,1)), length(ROIs_manual(:,1)) + day_i);
        plot(1:length(snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawPost{clustNum}(flight_i,1:xlimFlightPost)),snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawPost{clustNum}(flight_i,1:xlimFlightPost));
        hold on;
        if day_i == 1
            ylabel('Velocity (m/s)');
            xlabel('Time (s)');
            yt = get(gca,'YTick');
            xt = get(gca,'XTick');
            set(gca,'xticklabel',[]);
            %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        ylim([0 4.5]);
                title([num2str(length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})) ' flights']);
    end
    
    %plot all ROIs according to preferred sorting individually
    subplot(3,length(ROIs_manual(:,1)), 2*length(ROIs_manual(:,1)) + day_i);
    imagesc(snakeTraceT(day_i).snakeTrace.cRaw.normTracePost{clustNum}(snakeTraceT(day_i).snakeTrace.cRaw.IPost{clustNum},:),[1 4.5]);%1:xlimCalciumPost),[1 4.5]);
    colormap(hot);
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    else
        set(gca,'xticklabel',[]);
    end 
    hold off
end

% plot each ROI across all 9 days in same image with different images for each ROI sorted by day 1 preference
plotSnake_acrossDays_pref1_post = figure('units','normalized','outerposition',[0 0 0.5 1]);
sgtitle('Mean Post-Activity Across 9 Days for Same Flight Path (Pref Day 1): Gal 203011 to 200320');
col = jet(length(ROIs_manual(:,1)));
for day_i = 1:length(ROIs_manual(:,1))
    %plot flightpath for cluster #2
    for flight_i = 1:length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})
        subplot(3,length(ROIs_manual(:,1)),day_i);
        plot3(flightPathsF(day_i).flightPaths.pos(1,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(2,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(3,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),...
            '-','LineWidth',1,'Color', col(day_i,:));
        hold on;
        view(0,90);
        xlim([-3 3])
        ylim([-3 3])
        title(['Day ' num2str(day_i)]);
        if day_i == 1
            xlabel('m'); ylabel('m');
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        
        %plot velocity for cluster #2
        subplot(3,length(ROIs_manual(:,1)), length(ROIs_manual(:,1)) + day_i);
        plot(1:length(snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawPost{clustNum}(flight_i,1:xlimFlightPost)),snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawPost{clustNum}(flight_i,1:xlimFlightPost));
        hold on;
        if day_i == 1
            ylabel('Velocity (m/s)');
            xlabel('Time (s)');
            yt = get(gca,'YTick');
            xt = get(gca,'XTick');
            set(gca,'xticklabel',[]);
            %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        ylim([0 4.5]);
                title([num2str(length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})) ' flights']);
    end
    
    %plot all ROIs according to preferred sorting individually
    subplot(3,length(ROIs_manual(:,1)), 2*length(ROIs_manual(:,1)) + day_i);
    imagesc(snakeTraceT(day_i).snakeTrace.cRaw.normTracePost{clustNum}(snakeTraceT(1).snakeTrace.cRaw.IPost{clustNum},:),[1 4.5]); %1:xlimCalciumPost),[1 4.5]);
    colormap(hot);
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    else
        set(gca,'xticklabel',[]);
    end 
    hold off
end

%% plot all stable ROIs across all 9 days sorted by order of ROIs 
plotSnake_acrossDays_oddEven = figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle('Odd vs Even Trials Across 9 Days for Same Flight Path: Gal 203011 to 200320');
col = jet(length(ROIs_manual(:,1)));
for day_i = 1:length(ROIs_manual(:,1))
    %plot odd flightpath for cluster #2
    for flight_i = 1:2:length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})
        subplot(3,length(ROIs_manual(:,1))*2,(day_i*2)-1);
        plot3(flightPathsF(day_i).flightPaths.pos(1,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(2,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(3,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),...
            '-','LineWidth',1,'Color', col(day_i,:));
        hold on;
        view(0,90);
        xlim([-3 3])
        ylim([-3 3])
        title(['Day ' num2str(day_i)]);
        if day_i == 1
            xlabel('m'); ylabel('m');
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end   
        %plot odd velocity for cluster #2
        subplot(3,length(ROIs_manual(:,1))*2, length(ROIs_manual(:,1))*2 + (day_i*2)-1);
        plot(1:length(snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawFlight{clustNum}(flight_i,1:xlimFlight)),snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawFlight{clustNum}(flight_i,1:xlimFlight));
        hold on;
        if day_i == 1
            ylabel('Velocity (m/s)');
            xlabel('Time (s)');
            yt = get(gca,'YTick');
            xt = get(gca,'XTick');
            set(gca,'xticklabel',[]);
            %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        ylim([0 4.5]);
                title([num2str(length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})) ' flights']);
    end 
    %plot odd trials for all ROIs normalized to itself and plotted in 1-n numerical order
    %across all days
    subplot(3,2*length(ROIs_manual(:,1)), 4*length(ROIs_manual(:,1)) + (day_i*2)-1);
    imagesc(snakeTraceT(day_i).snakeTrace.cRaw.normTraceOdd{clustNum}(:,1:xlimCalcium),[1.5 4.5]);
    colormap(hot);
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    else
        set(gca,'xticklabel',[]);
    end 
    hold off
    
    %plot even flightpath for cluster #2
    for flight_i = 2:2:length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})
        subplot(3,length(ROIs_manual(:,1))*2,day_i*2);
        plot3(flightPathsF(day_i).flightPaths.pos(1,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(2,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),flightPathsF(day_i).flightPaths.pos(3,:,flightPathsF(day_i).flightPaths.clusterIndex{clustNum}(flight_i)),...
            '-','LineWidth',1,'Color', col(day_i,:));
        hold on;
        view(0,90);
        xlim([-3 3])
        ylim([-3 3])
        title(['Day ' num2str(day_i)]);
        if day_i == 1
            xlabel('m'); ylabel('m');
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end   
        %plot odd velocity for cluster #2
        subplot(3,length(ROIs_manual(:,1))*2, length(ROIs_manual(:,1))*2 + (day_i*2));
        plot(1:length(snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawFlight{clustNum}(flight_i,1:xlimFlight)),snakeTraceT(day_i).snakeTrace.cRaw.smoothSpeedRawFlight{clustNum}(flight_i,1:xlimFlight));
        hold on;
        if day_i == 1
            ylabel('Velocity (m/s)');
            xlabel('Time (s)');
            yt = get(gca,'YTick');
            xt = get(gca,'XTick');
            set(gca,'xticklabel',[]);
            %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
        else
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        ylim([0 4.5]);
                title([num2str(length(flightPathsF(day_i).flightPaths.clusterIndex{clustNum})) ' flights']);
    end 
    %plot odd trials for all ROIs normalized to itself and plotted in 1-n numerical order
    %across all days
    subplot(3,length(ROIs_manual(:,1))*2, length(ROIs_manual(:,1))*4 + (day_i*2));
    imagesc(snakeTraceT(day_i).snakeTrace.cRaw.normTraceEven{clustNum}(:,1:xlimCalcium),[1.5 4.5]);
    colormap(hot);
    if day_i == 1
        ylabel('ROI #');
        xlabel('Time (s)');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    else
        set(gca,'xticklabel',[]);
    end 
    hold off
end

%%save
if saveFlag ==1
   saveas(plotSnake_acrossDays,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays.tif']);
   savefig(plotSnake_acrossDays,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays.fig']);
   saveas(plotSnake_acrossDays_prefEach,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_prefEach.tif']);
   savefig(plotSnake_acrossDays_prefEach,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_prefEach.fig']);
   saveas(plotSnake_acrossDays_pref1,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_pref1.tif']);
   savefig(plotSnake_acrossDays_pref1,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_pref1.fig']);
   saveas(plotSnake_acrossDays_pre,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_pre.tif']);
   savefig(plotSnake_acrossDays_pre,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_pre.fig']);
   saveas(plotSnake_acrossDays_prefEach_pre,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_prefEach_pre.tif']);
   savefig(plotSnake_acrossDays_prefEach_pre,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_prefEach_pre.fig']);
   saveas(plotSnake_acrossDays_pref1_pre,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_pref1_pre.tif']);
   savefig(plotSnake_acrossDays_pref1_pre,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_pref1_pre.fig']);
   saveas(plotSnake_acrossDays_post,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_pre.tif']);
   savefig(plotSnake_acrossDays_post,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_pre.fig']);
   saveas(plotSnake_acrossDays_prefEach_post,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_prefEach_post.tif']);
   savefig(plotSnake_acrossDays_prefEach_post,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_prefEach_post.fig']);
   saveas(plotSnake_acrossDays_pref1_post,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_pref1_post.tif']);
   savefig(plotSnake_acrossDays_pref1_post,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_pref1_post.fig']);
   saveas(plotSnake_acrossDays_oddEven,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_oddEven.tif']);
   savefig(plotSnake_acrossDays_oddEven,[saveDir filesep 'Gal_200311to20_snakePlot_acrossDays_oddEven.fig']);
end