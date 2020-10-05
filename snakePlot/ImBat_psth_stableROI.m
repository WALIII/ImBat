function [psth_stableROI] = ImBat_psth_stableROI(varargin)
%function to plot firing fields as red dots against the flight paths of the
%bats for each day focusing only on the stable neurons from ROIs_manual
plotPSTHclustFlag = 0;
batId = '';
saveFlag = 0; %do you want to save the figures and output structure?
cRaw_flag = 1;
clustNum = 3; %which cluster trajectory to look at for each day

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'saveflag'
            saveFlag=varargin{i+1};
        case 'batid'
            batId = varargin{i+1};
        case 'plotpsthclustflag'
            plotPSTHclustFlag = varargin{i+1};
    end
end

if saveFlag == 1
saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
% Check if folder exists
if exist([saveDir1 datestr(now,'yymmdd') filesep 'psth_stableROI'])>0;
    disp('Youve been working today..');
else
    mkdir([saveDir1 datestr(now,'yymmdd') filesep 'psth_stableROI'])
end
saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'psth_stableROI' '\'];
end

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
    8 17 5 1 2 6 21 10 18 31 NaN 11 51 53 28 4 38 19 23 20;
    50 54 12 3 48 18 27 15 31 34 NaN NaN 28 NaN 29 25 24 22 38 14;
    8 NaN 4 28 3 18 10 35 42 25 13 NaN 50 39 46 NaN 49 2 32 26;
    14 NaN 3 28 2 6 33 26 18 45 NaN NaN 25 NaN 32 NaN 37 8 28 11];
g = dir('Ge*');
end
z = dir('Z1*');
dirTop = vertcat(g,z); %find all folders in top quality directory

%ROI_duplicate = cell(length(dirTop),1); %make cell for indices of duplicated ROIS

%initialize cells to collect vectors from each day
traceMean = cell(length(dirTop));
traceStd = cell(length(dirTop));
velocity = cell(length(dirTop));
flight = cell(length(dirTop));
traceIndiv = cell(length(dirTop));
colDay = jet(length(dirTop));

for day_i = 1:length(dirTop)
    %load results data
    try %extract metadata names and enter processed folder
        cd([dirTop(day_i).name filesep 'extracted'])
        flyFolders = dir('*fly*extraction');
        batName{day_i} = flyFolders(end).name(1:3);
        dateSesh{day_i} = flyFolders(end).name(5:10);
        sessionType{day_i} = flyFolders(end).name(12:16);
        cd(flyFolders(end).name);
    catch
        cd(dirTop(day_i).name);
        flyFolders = dir('*fly*extraction');
        batName{day_i} = flyFolders(end).name(1:3);
        dateSesh{day_i} = flyFolders(end).name(5:10);
        sessionType{day_i} = flyFolders(end).name(12:16);
        cd(flyFolders(end).name);
    end
    
    %load flightPaths and snakeTrace data
    %extract metadata names and enter analysis folder
    dirAnalysis = dir('analysis_*');
    if strcmp(batName{day_i}(1),'G')
        cd(dirAnalysis(end).name);
    else
        cd(dirAnalysis(end).name);
    end
    %load snakeTrace_cRaw data
    st = dir('*snakePlotData.mat');
    load(st(end).name);
    if cRaw_flag == 1
    snakeTrace_data = snakeTrace.cRaw;
    else
    snakeTrace_data = snakeTrace.s;
    end
    fp = dir('*flightPaths.mat');
    load(fp(end).name);
    close all;
    %take out means and stds from cluster 2 from each day and cell and save into variable
    %to plot later
    nROIs = size(ROIs_manual,2); %number of tractable ROIs
    nClusts = size(snakeTrace_data.smoothTraceRawFlight,2); %number of clusters in this day
    traceMean{day_i} = cell(nClusts,1);
    traceStd{day_i} = cell(nClusts,1);
    velocity{day_i} = cell(nClusts,1);
    flight{day_i} = cell(nClusts,1);
    traceIndiv{day_i} = cell(nClusts,1);
    
    for clust_ii = 1:nClusts
        nFlights = size(snakeTrace_data.smoothSpeedRawFlight{clust_ii},1);%size(flightPaths.clusterIndex{clust_ii},1);
        traceMean{day_i}{clust_ii} = zeros(nROIs,length(snakeTrace_data.normTraceFlight{clust_ii}(1,:)));
        traceStd{day_i}{clust_ii} = zeros(nROIs,length(snakeTrace_data.sdTraceFlight{clust_ii}(1,:)));
        velocity{day_i}{clust_ii} = cell(nROIs,1);
        flight{day_i}{clust_ii} = cell(nROIs,1);
        traceIndiv{day_i}{clust_ii} = cell(nROIs,1);
        
        for ROI_i = 1:nROIs
            try
                traceMean{day_i}{clust_ii}(ROI_i,:) = snakeTrace_data.normTraceFlight{clust_ii}(ROIs_manual(day_i,ROI_i),:);
                traceStd{day_i}{clust_ii}(ROI_i,:) = snakeTrace_data.sdTraceFlight{clust_ii}(ROIs_manual(day_i,ROI_i),:);
            catch
                traceMean{day_i}{clust_ii}(ROI_i,:) = zeros(1,length(snakeTrace_data.normTraceFlight{clust_ii}(1,:)));
                traceStd{day_i}{clust_ii}(ROI_i,:) = zeros(1,length(snakeTrace_data.sdTraceFlight{clust_ii}(1,:)));
            end
            velocity{day_i}{clust_ii}{ROI_i} = zeros(nFlights,length(snakeTrace_data.smoothSpeedRawFlight{clust_ii}(1,:)));
            flight{day_i}{clust_ii}{ROI_i} = zeros(nFlights,3,length(flightPaths.pos(1,:,1)));%,length(snakeTrace_cRaw.sdTraceFlight{clust_ii}(1,:)),length(snakeTrace_cRaw.sdTraceFlight{clust_ii}(1,:)));
            traceIndiv{day_i}{clust_ii}{ROI_i} = zeros(nFlights,length(snakeTrace_data.sdTraceFlight{clust_ii}(1,:)));
            for flight_i = 1:nFlights
                velocity{day_i}{clust_ii}{ROI_i}(flight_i,:) = snakeTrace_data.smoothSpeedRawFlight{clust_ii}(flight_i,:);
                flight{day_i}{clust_ii}{ROI_i}(flight_i,:,:,:) = flightPaths.pos(1:3,:,flightPaths.clusterIndex{clust_ii}(flight_i));
                try
                    traceIndiv{day_i}{clust_ii}{ROI_i}(flight_i,:) = snakeTrace_data.normTraceRawFlight{clust_ii}(flight_i,:,ROIs_manual(day_i,ROI_i));
                catch
                    traceIndiv{day_i}{clust_ii}{ROI_i}(flight_i,:) = zeros(1,length(snakeTrace_data.normTraceRawFlight{clust_ii}(flight_i,:,ROIs_manual(day_i,end))));
                end
            end
        end
    end
    
    %% plot the full psth with all clusters for tractable ROIs
    if plotPSTHclustFlag == 1
        %for each cell, for each cluster, plot flight paths, velocity, mean/std, and each trial
        nCells = size(ROIs_manual,2); %number of cells in this day
        nClusts = size(snakeTrace_data.smoothTraceRawFlight,2); %number of clusters in this day
        psth_each_cell = figure('units','normalized','outerposition',[0 0.1 0.8 0.8]);
        for cell_i = 1:nCells
            sgtitle(['Flight Aligned Activity per cluster: ' batName{day_i} ' ' dateSesh{day_i} ' Cell# ' num2str(cell_i)]);
            colClust = jet(nClusts); %color each cluster differently)
            try
                for clust_i = 1:nClusts
                    %plot flight paths of each cluster
                    nFlights = size(snakeTrace_data.smoothTraceRawFlight{clust_i},1); %number of flights in this cluster
                    colFlight = jet(nFlights); %color each cluster differently)
                    for flight_i = 1:nFlights
                        plotPaths = subplot(10,nClusts,[clust_i nClusts+clust_i]);
                        plot3(flightPaths.pos(1,:,flightPaths.clusterIndex{clust_i}(flight_i)),flightPaths.pos(2,:,flightPaths.clusterIndex{clust_i}(flight_i)),flightPaths.pos(3,:,flightPaths.clusterIndex{clust_i}(flight_i)),...
                            '-','LineWidth',1,'Color', colDay(clust_i,:));
                        hold on;
                        view(0,90);
                        xlim([-3 3])
                        ylim([-3 3])
                        title(['Clust ' num2str(clust_i)]);
                        if clust_i == 1
                            xlabel('m'); ylabel('m');
                        else
                            set(gca,'xticklabel',[],'yticklabel',[]);
                        end
                        
                        %plot velocity of each cluster
                        plotVelocity = subplot(10,nClusts,[clust_i+2*nClusts clust_i+3*nClusts]);
                        plot(1:length(snakeTrace_data.smoothSpeedRawFlight{clust_i}(flight_i,:)),snakeTrace_data.smoothSpeedRawFlight{clust_i}(flight_i,:),'Color',colDay(clust_i,:));
                        hold on;
                        if clust_i == 1
                            ylabel('Velocity (m/s)');
                            yt = get(gca,'YTick');
                            xt = get(gca,'XTick');
                            set(gca,'xticklabel',[]);
                            %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
                        else
                            set(gca,'xticklabel',[],'yticklabel',[]);
                        end
                        ylim([0 4.5]);
                        xlim([0 length(snakeTrace_data.smoothSpeedRawFlight{clust_i}(flight_i,:))]);
                        
                        %plot each flight activity
                        plotIndivTraces = subplot(10,nClusts,[clust_i+6*nClusts:nClusts:clust_i+9*nClusts]);
                        plot(1:length(snakeTrace_data.normTraceRawFlight{clust_i}(flight_i,:,ROIs_manual(day_i,cell_i))),(snakeTrace_data.normTraceRawFlight{clust_i}(flight_i,:,ROIs_manual(day_i,cell_i)))+flight_i*6,'Color',colDay(clust_i,:));
                        hold on;
                        if clust_i == 1
                            ylabel('Flight trials (df/f)');
                        end
                        xlabel('Time (s)');
                        xt = get(gca, 'XTick');
                        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
                        xlim([0 length(snakeTrace_data.normTraceRawFlight{clust_i}(flight_i,:,ROIs_manual(day_i,cell_i)))]);
                    end
                    %plot mean and stdev of activity per cluster
                    plotMeanTraces = subplot(10,nClusts,[clust_i+4*nClusts clust_i+5*nClusts]);
                    boundedline(1:length(psth_stableROI.traceMean{day_ii}{clustNum}(cell_ii,:)),psth_stableROI.traceMean{day_ii}{clustNum}(cell_ii,:),psth_stableROI.traceStd{day_ii}{clustNum}(cell_ii,:),...
                        'alpha','cmap',colDay(day_ii,:),'nan','gap');
                    hold on;
                    ylabel('Mean stdev df/f');
                    set(gca,'xticklabel',[]);
                    xlim([0 length(snakeTrace_data.normTraceFlight{clust_i}(ROIs_manual(day_i,cell_i),:))]);
                end
            catch
                sgtitle(['Flight Aligned Activity per cluster: ' batName{day_i} ' ' dateSesh{day_i} ' Cell# ' num2str(cell_i) ' NO ACTIVITY']);
            end
            
            % save
            % Save 'place cells' as jpg and fig files..
            if saveFlag == 1
                saveas(psth_each_cell,[saveDir filesep batName{day_i} '_' dateSesh{day_i} '_' sessionType{day_i} '_psth_clusts_' num2str(cell_i) '_' datestr(now,'yymmdd-HHMM') '.tif']);
                savefig(psth_each_cell,[saveDir filesep batName{day_i} '_' dateSesh{day_i} '_' sessionType{day_i} '_psth_clusts_' num2str(cell_i) '_' datestr(now,'yymmdd-HHMM') '.fig']);
                %saveas(psth_each_cell,[pwd filesep batName '_' dateSesh '_' sessionType '_psth_clusts' num2str(cell_i) '.svg']);
            end
            clf;
            
        end        
    end
    cd(dirTop(day_i).folder)
end

    psth_stableROI.traceMean = traceMean;
    psth_stableROI.traceStd = traceStd;
    psth_stableROI.velocity = velocity;
    psth_stableROI.flight = flight;
    psth_stableROI.traceIndiv = traceIndiv;
    if saveFlag ==1
        save([saveDir filesep batName{day_i} '_' dateSesh{day_i} '_' sessionType{day_i} 'psth_stableROIdata_' datestr(now,'yymmdd-HHMM') '.mat'],'psth_stableROI');
    end

%%
plotPSTH_acrossDays = figure('units','normalized','outerposition',[0 0.1 0.4 0.8]);
for cell_ii = 1:nROIs
    for day_ii = 1:length(dirTop)
        try
            %plot flight paths of each cluster
            nFlights = size(psth_stableROI.velocity{day_ii}{clustNum}{1},1); %number of flights in this cluster
            colFlight = jet(nFlights); %color each cluster differently)
            for flight_i = 1:nFlights
                plotPaths = subplot(10,1,[1 2]);
                plot3(squeeze(psth_stableROI.flight{day_ii}{clustNum}{cell_ii}(flight_i,1,:)),squeeze(psth_stableROI.flight{day_ii}{clustNum}{cell_ii}(flight_i,2,:)),squeeze(psth_stableROI.flight{day_ii}{clustNum}{cell_ii}(flight_i,3,:)),...
                    '-','LineWidth',1,'Color',colDay(day_ii,:));
                hold on;
                view(0,90);
                xlim([-3 3])
                ylim([-3 3])
                if strcmp(batId,'Gal')
                    title(['ROI # ' num2str(cell_ii) ': Gal 200311to20 - clust' num2str(clustNum)]);
                elseif strcmp(batId,'Gen')
                    title(['ROI # ' num2str(cell_ii) ': Gen 200319to24 - clust' num2str(clustNum)]);
                end
                xlabel('m'); ylabel('m');
                
                %plot velocity of each cluster
                plotVelocity = subplot(10,1,[3 4]);
                plot(1:length(psth_stableROI.velocity{day_ii}{clustNum}{cell_ii}(flight_i,:)),psth_stableROI.velocity{day_ii}{clustNum}{cell_ii}(flight_i,:),'Color',colDay(day_ii,:));
                hold on;
                ylabel('Velocity (m/s)');
                yt = get(gca,'YTick');
                xt = get(gca,'XTick');
                set(gca,'xticklabel',[]);
                %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
                ylim([0 4.5]);
                xlim([0 length(psth_stableROI.velocity{1}{clustNum}{cell_ii}(1,:))]);
                
            end
            %plot mean and stdev of activity per cluster
            plotMeanTraces = subplot(10,1,[5 10]);
            p(day_ii) = boundedline(1:length(psth_stableROI.traceMean{day_ii}{clustNum}(cell_ii,:)),psth_stableROI.traceMean{day_ii}{clustNum}(cell_ii,:),psth_stableROI.traceStd{day_ii}{clustNum}(cell_ii,:),...
                'alpha','cmap',colDay(day_ii,:),'nan','gap');
            hold on;
            ylabel('Mean stdev df/f');
            xlabel('Time (s)');
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlim([0 length(psth_stableROI.traceMean{1}{clustNum}(cell_ii,:))]);
        catch
            sgtitle(['Flight Aligned Activity per cluster: ' batName{day_ii} ' ' dateSesh{day_ii} ' Cell# ' num2str(cell_ii) ' NO ACTIVITY']);
        end
        
    end
    legend([p(1),p(2),p(3),p(4),p(5)],'Day 1','Day 2','Day 3','Day 4','Day 5');

    drawnow;
    if saveFlag == 1
        if strcmp(batId,'Gal')
        saveas(plotPSTH_acrossDays,[saveDir filesep 'Gal_200311to20_psth_acrossDays_clust' num2str(clustNum) '_ROI' num2str(cell_ii) '_' datestr(now,'yymmdd_HHMM') '.tif']);
        savefig(plotPSTH_acrossDays,[saveDir filesep 'Gal_200311to20_psth_acrossDays_clust' num2str(clustNum) '_ROI' num2str(cell_ii) '_' datestr(now,'yymmdd_HHMM') '.fig']);
        elseif strcmp(batId,'Gen')
        saveas(plotPSTH_acrossDays,[saveDir filesep 'Gen_200319to24_psth_acrossDays_clust' num2str(clustNum) '_ROI' num2str(cell_ii) '_' datestr(now,'yymmdd_HHMM') '.tif']);
        savefig(plotPSTH_acrossDays,[saveDir filesep 'Gen_200319to24_psth_acrossDays_clust' num2str(clustNum) '_ROI' num2str(cell_ii) '_' datestr(now,'yymmdd_HHMM') '.fig']);
        end
    end
    clf;
end
cd(dirTop(day_i).folder)
end
%
%
%
% % %% load all figures from 1 ROI and save in subplots so can compare across days for same ROI
% % cd(['plots' filesep datestr(now,'yymmdd')]);
% % for cell_ii = 1:length(ROIs_manual(1,:))
% %     figDir = dir(['*psth_clusts_' num2str(cell_ii) '.fig'])
% %     plotPSTHAcrossDays(cell_ii) = figure('units','normalized','outerposition',[0 0 1 0.8]);
% %     sgtitle(['Gal PSTH 200311-200320: Clust2 ROI' num2str(cell_ii)]);
% %     for ii = 1:length(figDir)
% %         figure(plotPSTHAcrossDays(cell_ii))
% %         % Prepare subplots
% %                 velocities(ii)=subplot(10,1,3:4);
% %         meanTraces(ii)=subplot(10,1,5:10);
% %         paths(ii)=subplot(10,1,1:2);
% %         title(['Day ' num2str(ii)]);
% %         hold on;
% %
% %         % Load saved figures
% %         figPaths(ii)=hgload(figDir(ii).name);
% %         % Paste figures on the subplots
% %         copyobj(figPaths(ii).plotPaths,paths(ii));
% %         copyobj(figPaths(ii).plotVelocity,velocities(ii));
% %         copyobj(figPaths(ii).plotMeanTraces,meanTraces(ii));
% %         %copyobj(allchild(get(figPaths(ii),'CurrentAxes')),paths(ii));
% %
% %     end
% %     if saveFlag == 1
% %         saveas(plotFiringTrajectoryAcrossDays(cell_ii),[pwd filesep datestr(now,'yymmdd-HHMM') '_' 'Gal_200311to20_placeCell_' num2str(cell_ii) '_acrossDays.tif']);
% %         savefig(plotFiringTrajectoryAcrossDays(cell_ii),[pwd filesep datestr(now,'yymmdd-HHMM') '_' 'Gal_200311to20_placeCell_' num2str(cell_ii) '_acrossDays.fig']);
% %     end
% %
% %     close gcf;
% % end

