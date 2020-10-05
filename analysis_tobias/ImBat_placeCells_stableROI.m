function ImBat_placeCells_stableROI
%function to plot firing fields as red dots against the flight paths of the
%bats for each day focusing only on the stable neurons from ROIs_manual

saveFlag = 1; %do you want to save the figures and output structure?
saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
% Check if folder exists
if exist([saveDir1 datestr(now,'yymmdd')])>0;
    disp('Youve been working today..');
else
    mkdir([saveDir1 datestr(now,'yymmdd')])
end
saveDir = [saveDir1 datestr(now,'yymmdd') '\'];

offset = 0; % account for slow calcium estimation ~move locations back 100ms in time... This is the knob to turn for 'prospective' coding...
speedThresh = 0.7; %threshold for when to eliminate nonflying moments
spikeThreshMult = 5; %number of times to multiply std of s vector for determining spike threshold
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

%ROI_duplicate = cell(length(dirTop),1); %make cell for indices of duplicated ROIS


for d = 1:length(dirTop)
    %load results data
    try %extract metadata names and enter processed folder
        cd([dirTop(d).name filesep 'extracted'])
        flyFolders = dir('*fly*extraction');
        batName{d} = flyFolders(end).name(1:3);
        dateSesh{d} = flyFolders(end).name(5:10);
        sessionType{d} = flyFolders(end).name(12:16);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        if strcmp(batName{d}(1),'G')
            cd(dirProcessed(end).name);
        else
            cd(dirProcessed(end).name);
        end
    catch
        cd(dirTop(d).name);
        flyFolders = dir('*fly*extraction');
        batName{d} = flyFolders(end).name(1:3);
        dateSesh{d} = flyFolders(end).name(5:10);
        sessionType{d} = flyFolders(end).name(12:16);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        cd(dirProcessed(end).name);
    end
    cellData = load('results.mat');
    alignment = load('Alignment.mat');
    cd(dirProcessed(end).folder);
    
    %load flightPaths data
    %extract metadata names and enter analysis folder
    dirAnalysis = dir('analysis_*');
    if strcmp(batName{d}(1),'G')
        cd(dirAnalysis(end).name);
    else
        cd(dirAnalysis(end).name);
    end
    fp = dir('*flightPaths.mat');
    load(fp(end).name);
    close all;
    
    %%
    % Remove that pesky first flight that starts at 0:0;
    try
        a = find(alignment.out.flights(:,1) == 0);
        a2 = isnan(alignment.out.flights(a(1):end,1));
        a3 = find(a2>0);
        alignment.out.flights(a(1):a(1)+a3(1),:) = NaN;
        alignment.out.Location2 = alignment.out.flights;
    catch
    end
    
    % Plot the location in space that each cell is active in 1 figure
    plotFiringTrajectory =  figure('units','normalized','outerposition',[0 0 1 0.8]);
    sgtitle([batName{d} ' ' dateSesh{d} ' ' sessionType{d} ': Firing Fields']);
    n = 1;
    for ii = ROIs_manual(d,:)%1:length(cellData.results.S(:,1)); % for each cell
        set(0,'CurrentFigure',plotFiringTrajectory);
        subplot(ceil(length(ROIs_manual(d,:))/3),3,n)
        hold on;
        plot(alignment.out.flights(:,1),alignment.out.flights(:,2),'k');% plot the flight trajectory in space
        %plot3(alignment.out.flights(:,1),alignment.out.flights(:,2),alignment.out.flights(:,3),'k');%,'LineWidth',2);% plot the flight trajectory in space
        try
            spikeThresh(n) = median(cellData.results.S(ii,:)) + std(cellData.results.S(ii,:))*spikeThreshMult; %threshold for eliminating noise from the S matrix
            [~,xy] = find(cellData.results.S(ii,:)>spikeThresh(n));  % get time neuron is active
            xy = xy(xy<length(alignment.out.video_times));
            Spike_times = alignment.out.video_times(xy)-offset; % convert this to 'spike time'
            
            %only take data from when bat is flying
            flightVect = alignment.out.flights; %make a new flight vector to eliminate the subthreshold flight_times
            [yz,~] = find(flightPaths.batSpeed<speedThresh); %find when bat is not flying
            flightVect(yz) = NaN; %set nonflying times to NaN
            peak_heights = cellData.results.S(ii,xy);
            
            LX = zeros(1,length(Spike_times(:,1)));
            LY = zeros(1,length(Spike_times(:,1)));
            LZ = zeros(1,length(Spike_times(:,1)));
            PH = zeros(1,length(Spike_times(:,1)));
            
            try % this 'try/catch' is to avoid crashing if cells are not active in plotting window...
                for i = 1:size(Spike_times,1)
                    try
                        % Find the closest 'Location time' to the 'Spike time'
                        [minValue(:,i),closestIndex(:,i)] = min(abs(alignment.out.Location_time-Spike_times(i)));
                        LX(i) = flightVect(closestIndex(:,i),1);%alignment.out.flights(closestIndex(:,i),1);
                        LY(i) = flightVect(closestIndex(:,i),2);%alignment.out.flights(closestIndex(:,i),2);
                        LZ(i) = flightVect(closestIndex(:,i),3);%alignment.out.flights(closestIndex(:,i),3);
                        if isnan(LX(i))
                            PH(i) = NaN; %make the peak into a NaN for scaling purposes
                        else
                            PH(i) = full(peak_heights(i));
                        end
                    catch % we need this if Spiketime occurs before/after the location tracking was on..
                        disp('cell catch');
                        continue
                    end
                end
                % display, in the title, how many bursts there were:
                LX_s = LX;
                LX_s(isnan(LX)) = [];
                disp([num2str(size(LX_s)),' Bursts in flight'])
                PH = mat2gray(PH);
                hold on
                %scatter(LX,LY,(PH*75)+1,'or','filled');
                scatter(LX,LY,(PH*75)+1,'or','filled');
                %scatter3(LX,LY,LZ,(PH*75)+1,'or','filled');
                %uistack(dots,'top');
                % % modify labels for tick marks
                xlim([-3000 3000]);
                ylim([-3000 3000]);
                set(gca,'xticklabel',[],'yticklabel',[]);
                title(['ROI ' num2str(n) '(' num2str(ii)  '): ' num2str(size(LX_s)) ' Bursts']);
                %hold off
            catch % if cell was not active...
                disp('cell not active');
                continue
            end
            
        catch
            xlim([-3000 3000]);
            ylim([-3000 3000]);
            set(gca,'xticklabel',[],'yticklabel',[]);
            title(['ROI ' num2str(n) '(' num2str(ii)  '): Not Active']);
        end
        
        % Plot the location in space that each cell is active in individual figures
        plotFiringTrajectoryIndiv =  figure();
        a2 = axes;
        try
            plot(a2,alignment.out.flights(:,1),alignment.out.flights(:,2),'k');% plot the flight trajectory in space
            hold on;
            disp([num2str(size(LX_s)),' Bursts in flight'])
            PH = mat2gray(PH);
            %scatter(LX,LY,(PH*75)+1,'or','filled');
            scatter(a2,LX,LY,(PH*75)+1,'or','filled');
            %scatter3(LX,LY,LZ,(PH*75)+1,'or','filled');
            %uistack(dots,'top');
            % % modify labels for tick marks
            xlim([-3000 3000]);
            ylim([-3000 3000]);
            xticks = get(gca,'xtick');
            yticks = get(gca,'ytick');
            newlabelsX = arrayfun(@(ax) sprintf('%g', ax/1000), xticks, 'un', 0);
            newlabelsY = arrayfun(@(ay) sprintf('%g', ay/1000), yticks, 'un', 0);
            set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
            xlabel('m'); ylabel('m');
            title([batName{d} ' ' dateSesh{d} ' ' sessionType{d} ':ROI ' num2str(n) '(' num2str(ii) '): ' num2str(size(LX_s)) ' Bursts']);
            hold off
            
        catch % if cell was not active...
            disp('cell not active');
            xlim([-3000 3000]);
            ylim([-3000 3000]);
            xticks = get(gca,'xtick');
            yticks = get(gca,'ytick');
            newlabelsX = arrayfun(@(ax) sprintf('%g', ax/1000), xticks, 'un', 0);
            newlabelsY = arrayfun(@(ay) sprintf('%g', ay/1000), yticks, 'un', 0);
            set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
            xlabel('m'); ylabel('m');
            title([batName{d} ' ' dateSesh{d} ' ' sessionType{d} ':ROI ' num2str(n) '(' num2str(ii) '): Not Active']);
            hold off
        end
        % Save 'place cells' as jpg and fig files..
        if saveFlag == 1
            %set(findall(gcf,'-property','FontSize'),'FontSize',20);
            saveas(plotFiringTrajectoryIndiv,[saveDir filesep datestr(now,'yymmdd-HHMM') '_' batName{d} '_' dateSesh{d} '_' sessionType{d} '_placeCell_' num2str(n) '.tif']);
            savefig(plotFiringTrajectoryIndiv,[saveDir filesep datestr(now,'yymmdd-HHMM') '_' batName{d} '_' dateSesh{d} '_' sessionType{d} '_placeCell_' num2str(n) '.fig']);
            %saveas(gcf,[saveDir '\' batName '_' dateSesh '_' sessionType '_placeCell_' num2str(n) '.svg']);
        end
        
        % Clear the buffer for the next cell:
        clear LX LY LZ closestIndex Spike_times
        n = n +1;
        
        %clf
    end
    
    % Save 'place cells' as jpg and fig files..
    if saveFlag == 1
        %set(findall(gcf,'-property','FontSize'),'FontSize',20);
        saveas(plotFiringTrajectory,[saveDir '\' datestr(now,'yymmdd-HHMM') '_' batName{d} '_' dateSesh{d} '_' sessionType{d} '_placeCell_all.tif']);
        savefig(plotFiringTrajectory,[saveDir filesep datestr(now,'yymmdd-HHMM') '_' batName{d} '_' dateSesh{d} '_' sessionType{d} '_placeCell_all.fig']);
        %saveas(gcf,[saveDir '\' batName '_' dateSesh '_' sessionType '_placeCell_all.svg']);
    end
    
    
    cd(dirTop(d).folder)
end
close all;
%% load all figures from 1 ROI and save in subplots so can compare across days for same ROI

cd(['plots' filesep datestr(now,'yymmdd')]);
for i = 1:length(ROIs_manual(1,:))
    figDir = dir(['*placeCell_' num2str(i) '.fig'])
    plotFiringTrajectoryAcrossDays(i) = figure('units','normalized','outerposition',[0 0 1 0.8]);
    sgtitle(['Gal ROI Selectivity 200311-200320: ROI ' num2str(i)]);
    for ii = 1:length(figDir)
        figure(plotFiringTrajectoryAcrossDays(i))
        h(ii)=subplot(ceil(length(ROIs_manual(:,ii))/3),3,ii);
        title(['Day ' num2str(ii)]);
        hold on;
        % Load saved figures
        figL(ii)=hgload(figDir(ii).name);
        % Prepare subplots
        % Paste figures on the subplots
        copyobj(allchild(get(figL(ii),'CurrentAxes')),h(ii));
        
    end
    if saveFlag == 1
    saveas(plotFiringTrajectoryAcrossDays(i),[pwd filesep datestr(now,'yymmdd-HHMM') '_' 'Gal_200311to20_placeCell_' num2str(i) '_acrossDays.tif']);
    savefig(plotFiringTrajectoryAcrossDays(i),[pwd filesep datestr(now,'yymmdd-HHMM') '_' 'Gal_200311to20_placeCell_' num2str(i) '_acrossDays.fig']);
    end
    
    close all;
end

cd(dirTop(d).folder)

