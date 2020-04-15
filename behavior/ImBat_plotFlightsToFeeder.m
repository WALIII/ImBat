function [flightPathsToFeeder, flightPathsFromFeeder, flightPathsClusterToFeederEach, flightPathsClusterToFeederAll, flightFeedersStartStop]...
    = ImBat_plotFlightsToFeeder(trackData,varargin)

%length of time for flight up until reaching feeder (s)
preflightTime = 4;
%length of time to wait before bat leaves feeder for FlightFromFeeder track %(s)
delay = 3;
%number of clusters for kmeans clustering of flight trajectories
nclusters = 6;
ntrajectories = 6; %number of output trajectories from kmeans that you want to look at

batName = [];
dateSesh = [];
sessionType = [];
saveFlag = 0;

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'batname'
            batName=varargin{i+1};
        case 'datesesh'
            dateSesh = varargin{i+1};
        case 'sessiontype'
            sessionType = varargin{i+1};
        case 'saveflag'
            saveFlag = varargin{i+1};
        case 'analysisfolder'
            analysis_Folder = varargin{i+1};
    end
end

%labels for loading and saving data if running independent fromImBat_analyze
if saveFlag == 1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
    
    cellData = load([pwd '/processed/Motion_corrected_Data_DS_results.mat']);
    alignment = load([pwd '/processed/Alignment.mat']);
    load([pwd '/' analysis_Folder '/' label '_flightPaths.mat']);
end
 



%find events when bat landed on feeder
event_ttls = trackData.AnalogSignals(:,1); %trial data
[R,LT,UT,LL,UL] = risetime(event_ttls,trackData.VideoFrameRate);
%markerSet = trackData.Markers;

mx = trackData.Markers(:,1,1);
my = trackData.Markers(:,1,2);
mz = trackData.Markers(:,1,3);

%set zeros to nan
mx(find(mx == 0)) = nan;
my(find(my == 0)) = nan;
mz(find(mz == 0)) = nan;

mxFull = fillmissing(mx,'nearest');
myFull = fillmissing(my,'nearest');
mzFull = fillmissing(mz,'nearest');

%threshold based on speed
Vx = gradient(mxFull, 1/trackData.VideoFrameRate);
Vy = gradient(myFull, 1/trackData.VideoFrameRate);
Vz = gradient(mzFull, 1/trackData.VideoFrameRate);
batSpeed = sqrt(Vx.^2 + Vy.^2 + Vz.^2)/1000; %in m/s
batSpeedMean = movmean(batSpeed,8);
batSpeedFilt = medfilt1(batSpeedMean,20);

%batSpeedMean(find(batSpeedMean>4.5)) = nan;

%% plot flight paths of bat to feeder
flightPathsToFeeder = figure();
CM = jet(length(LT));
flightToFeederEnd = zeros(length(LT),1);
flightToFeederStart = zeros(length(LT),1);
%try
for i = 1:length(LT)
    flightToFeederEnd(i) = round(LT(i) * trackData.VideoFrameRate);
    flightToFeederStart(i) = find(batSpeedMean(1:flightToFeederEnd(i))==0,1,'last'); %minus the preflight # of seconds flightToFeederEnd(i) - preFlightTime
    %scatter3(markerSet(flightStart:flightEnd,1,1),markerSet(flightStart:flightEnd,1,2),markerSet(flightStart:flightEnd,1,3))
    plot3(mxFull(flightToFeederStart(i):flightToFeederEnd(i)),myFull(flightToFeederStart(i):flightToFeederEnd(i)),mzFull(flightToFeederStart(i):flightToFeederEnd(i)),'LineWidth',2,'color',CM(i,:))
    hold on
    %get start/stop xyz position of the flights to feeders
    fFeedStartxyz(i,1) = round(nanmean(mxFull(flightToFeederStart(i):flightToFeederStart(i))));%+trackData.VideoFrameRate/4)));
    fFeedStartxyz(i,2) = round(nanmean(myFull(flightToFeederStart(i):flightToFeederStart(i))));%+trackData.VideoFrameRate/4)));
    fFeedStartxyz(i,3) = round(nanmean(mzFull(flightToFeederStart(i):flightToFeederStart(i))));%+trackData.VideoFrameRate/4)));
    
    fFeedEndxyz(i,1) = round(nanmean(mxFull(flightToFeederEnd(i):flightToFeederEnd(i)+trackData.VideoFrameRate/2)));
    fFeedEndxyz(i,2) = round(nanmean(myFull(flightToFeederEnd(i):flightToFeederEnd(i)+trackData.VideoFrameRate/2)));
    fFeedEndxyz(i,3) = round(nanmean(mzFull(flightToFeederEnd(i):flightToFeederEnd(i)+trackData.VideoFrameRate/2)));
    
    scatter3(fFeedStartxyz(i,1),fFeedStartxyz(i,2),fFeedStartxyz(i,3),50,'r','filled')
    hold on
    scatter3(fFeedEndxyz(i,1),fFeedEndxyz(i,2),fFeedEndxyz(i,3),50,'k','filled')
end
title(['Flights to feeder: ' batName ' ' dateSesh ' ' sessionType]);
% % modify labels for tick marks
view(0,90);
xlim([-3000 3000]);
ylim([-3000 3000]);
xticks = get(gca,'xtick');
yticks = get(gca,'ytick');
newlabelsX = arrayfun(@(ax) sprintf('%g', ax/1000), xticks, 'un', 0);
newlabelsY = arrayfun(@(ay) sprintf('%g', ay/1000), yticks, 'un', 0);
set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
xlabel('m'); ylabel('m');
hold off

%% plot flight paths of bat from feeder
flightFromFeederEnd = zeros(length(LT));
flightFromFeederStart = zeros(length(LT));
flightPathsFromFeeder = figure();
for i = 1:length(LT)
    flightFromFeederMid(i) = find(batSpeedFilt(flightToFeederEnd(i):end)>3,1); %minus the preflight # of seconds flightToFeederEnd(i) - preFlightTime
    %flightFromFeederStart(i) = round(LT(i) * trackData.VideoFrameRate);
    flightFromFeederEnd(i) = flightToFeederEnd(i) + flightFromFeederMid(i) + find(batSpeedMean(flightToFeederEnd(i)+flightFromFeederMid(i):end)==0,1);
    %flightFromFeederEnd(i) = flightFromFeederStart(i) + delay + (trackData.VideoFrameRate*preflightTime);
    %check that the end of the flight does not occur after session finishes
    if flightFromFeederEnd(i) > length(mxFull)
        flightFromFeederEnd(i) = length(mxFull)
    end
    plot3(mxFull(flightToFeederEnd(i):flightFromFeederEnd(i)),myFull(flightToFeederEnd(i):flightFromFeederEnd(i)),mzFull(flightToFeederEnd(i):flightFromFeederEnd(i)),'LineWidth',2,'color',CM(i,:))
    hold on
        %get start/stop xyz position of the flights to feeders
    %fLeaveStartxyz(i,1) = round(nanmean(mxFull(flightFromFeederStart(i):flightFromFeederStart(i)+trackData.VideoFrameRate/4)));
    %fLeaveStartxyz(i,2) = round(nanmean(myFull(flightFromFeederStart(i):flightFromFeederStart(i)+trackData.VideoFrameRate/4)));
    %fLeaveStartxyz(i,3) = round(nanmean(mzFull(flightFromFeederStart(i):flightFromFeederStart(i)+trackData.VideoFrameRate/4)));
    
    fLeaveEndxyz(i,1) = round(nanmean(mxFull(flightFromFeederEnd(i):flightFromFeederEnd(i))));%+trackData.VideoFrameRate/2)));
    fLeaveEndxyz(i,2) = round(nanmean(myFull(flightFromFeederEnd(i):flightFromFeederEnd(i))));%+trackData.VideoFrameRate/2)));
    fLeaveEndxyz(i,3) = round(nanmean(mzFull(flightFromFeederEnd(i):flightFromFeederEnd(i))));%+trackData.VideoFrameRate/2)));
    
    scatter3(fFeedEndxyz(i,1),fFeedEndxyz(i,2),fFeedEndxyz(i,3),50,'r','filled')
    hold on
    scatter3(fLeaveEndxyz(i,1),fLeaveEndxyz(i,2),fLeaveEndxyz(i,3),50,'k','filled')
    %pause(0.1)
end
title(['Flights from feeder: ' batName ' ' dateSesh ' ' sessionType]);
 % % modify labels for tick marks
        view(0,90); 
        xlim([-3000 3000]);
        ylim([-3000 3000]);        
        xticks = get(gca,'xtick');
        yticks = get(gca,'ytick');
        newlabelsX = arrayfun(@(ax) sprintf('%g', ax/1000), xticks, 'un', 0);
        newlabelsY = arrayfun(@(ay) sprintf('%g', ay/1000), yticks, 'un', 0);
        set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
        xlabel('m'); ylabel('m');
        hold off

%% k means cluster of flight trajectories into nclusters
%find pairs of start and endpoints with a high number of flights

rng(2) %control random number generation
try
kstart = kmeans(fFeedStartxyz,nclusters);
rng(2)
kend = kmeans(fFeedEndxyz,nclusters);
catch
kstart = kmeans(fFeedStartxyz,nclusters/2);
rng(2)
kend = kmeans(fFeedEndxyz,nclusters/2);
end
nflights = [];
allflights = [];
npair = [];
ncounter = 0;
for ks = 1 : nclusters
    for ke = 1 : nclusters
        ncounter = ncounter + 1;
        npair(ncounter,:) = [ks ke];
        nflights(ncounter) = size(intersect(find(kstart == ks),find(kend == ke)),1);
        allflights{ncounter} = intersect(find(kstart == ks),find(kend == ke))';
    end
end
 
[~, ssf] = sort(nflights,'descend'); %sort clustered flights

%% plot clustered flights individually by each cluster
flightPathsClusterToFeederEach = figure();
jj = jet(ntrajectories);
for traj = 1 : 10
    try
    if traj<ntrajectories + 1 % Only plot the first ntrajectories...
    subplot(round(ntrajectories/2),2,traj)
    for nf = allflights{ssf(traj)}
        plot3(mxFull(flightToFeederStart(nf):flightToFeederEnd(nf)),myFull(flightToFeederStart(nf):flightToFeederEnd(nf)),mzFull(flightToFeederStart(nf):flightToFeederEnd(nf)),'LineWidth',1,'Color',jj(traj,:))
        hold on
        scatter3(fFeedStartxyz(nf,1),fFeedStartxyz(nf,2),fFeedStartxyz(nf,3),50,'r','filled')
        %hold on
        scatter3(fFeedEndxyz(nf,1),fFeedEndxyz(nf,2),fFeedEndxyz(nf,3),50,'k','filled')
        hold on
 % % modify labels for tick marks
        view(0,90); 
        xlim([-3000 3000]);
        ylim([-3000 3000]);        
        xticks = get(gca,'xtick');
        yticks = get(gca,'ytick');
        newlabelsX = arrayfun(@(ax) sprintf('%g', ax/1000), xticks, 'un', 0);
        newlabelsY = arrayfun(@(ay) sprintf('%g', ay/1000), yticks, 'un', 0);
        set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
        xlabel('m'); ylabel('m');
        hold off
    end
    end
        clusterIndex{traj}(:) = allflights{ssf(traj)};

    catch
        disp(' no more flights...');
    end
end
sgtitle(['Flight Clusters to Feeder start(r)/stop(b): ' batName ' ' dateSesh ' ' sessionType]);
%xlabel('mm'); ylabel('mm'); zlabel('mm');
hold off

%% plot all the clusters in 1 figure in different colors
flightPathsClusterToFeederAll = figure();
jj = jet(ntrajectories);
for traj = 1 : ntrajectories
    for nf = allflights{ssf(traj)}
        plot3(mxFull(flightToFeederStart(nf):flightToFeederEnd(nf)),myFull(flightToFeederStart(nf):flightToFeederEnd(nf)),mzFull(flightToFeederStart(nf):flightToFeederEnd(nf)),'LineWidth',1,'Color',jj(traj,:))
        hold on
        scatter3(fFeedStartxyz(nf,1),fFeedStartxyz(nf,2),fFeedStartxyz(nf,3),50,'r','filled')
        hold on
        scatter3(fFeedEndxyz(nf,1),fFeedEndxyz(nf,2),fFeedEndxyz(nf,3),50,'k','filled')
        hold on
    end
end
title(['Flight Clusters to Feeder start(r)/stop(b): ' batName ' ' dateSesh ' ' sessionType]);
% % modify labels for tick marks
view(0,90);
xlim([-3000 3000]);
ylim([-3000 3000]);
xticks = get(gca,'xtick');
yticks = get(gca,'ytick');
newlabelsX = arrayfun(@(ax) sprintf('%g', ax/1000), xticks, 'un', 0);
newlabelsY = arrayfun(@(ay) sprintf('%g', ay/1000), yticks, 'un', 0);
set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
xlabel('m'); ylabel('m');
hold off

flightFeedersStartStop.clusterIndex = clusterIndex;
flightFeedersStartStop.flightToFeederEnd = flightToFeederEnd;
flightFeedersStartStop.flightToFeederStart = flightToFeederStart;
flightFeedersStartStop.flightFromFeederEnd = flightFromFeederEnd;
flightFeedersStartStop.flightFromFeederStart = flightFromFeederStart;
flightFeedersStartStop.fFeedStartxyz = fFeedStartxyz;
flightFeedersStartStop.fFeedEndxyz = fFeedEndxyz;
% catch
%     
%     flightFeedersStartStop.clusterIndex = [];
%     flightPathsFromFeeder = figure();
%     flightPathsClusterToFeederEach = figure();
%     flightPathsClusterToFeederAll = figure();
%     flightFeedersStartStop.flightToFeederEnd = [];
%     flightFeedersStartStop.flightToFeederStart = [];
%     flightFeedersStartStop.flightFromFeederEnd = [];
%     flightFeedersStartStop.flightFromFeederStart = [];
%     flightFeedersStartStop.fFeedStartxyz = [];
%     flightFeedersStartStop.fFeedEndxyz = [];
%     
% end

beta = 1;
