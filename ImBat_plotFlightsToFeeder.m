function [flightPathsToFeeder, flightPathsFromFeeder, flightPathsClusterToFeederEach, flightPathsClusterToFeederAll, flightFeedersStartStop]...
    = ImBat_plotFlightsToFeeder(trackData)

global batName dateSesh sessionType

%length of time for flight up until reaching feeder (s)
preflightTime = 2;
%length of time to wait before bat leaves feeder for FlightFromFeeder track %(s)
delay = 3;
%number of clusters for kmeans clustering of flight trajectories
nclusters = 4;
ntrajectories = 6; %number of output trajectories from kmeans that you want to look at 



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

mx = fillmissing(mx,'nearest');
my = fillmissing(my,'nearest');
mz = fillmissing(mz,'nearest');

%plot flight paths of bat to feeder
flightPathsToFeeder = figure();
CM = jet(length(LT));
flightToFeederEnd = zeros(length(LT));
flightToFeederStart = zeros(length(LT));
try
for i = 1:length(LT)
    flightToFeederEnd(i) = round(LT(i) * trackData.VideoFrameRate);
    flightToFeederStart(i) = flightToFeederEnd(i) - (trackData.VideoFrameRate*preflightTime); %minus the preflight # of seconds
    %scatter3(markerSet(flightStart:flightEnd,1,1),markerSet(flightStart:flightEnd,1,2),markerSet(flightStart:flightEnd,1,3))
    plot3(mx(flightToFeederStart(i):flightToFeederEnd(i)),my(flightToFeederStart(i):flightToFeederEnd(i)),mz(flightToFeederStart(i):flightToFeederEnd(i)),'LineWidth',2,'color',CM(i,:))
    hold on
    %get start/stop xyz position of the flights to feeders
    fFeedStartxyz(i,1) = round(nanmean(mx(flightToFeederStart(i):flightToFeederStart(i)+trackData.VideoFrameRate/2)));
    fFeedStartxyz(i,2) = round(nanmean(my(flightToFeederStart(i):flightToFeederStart(i)+trackData.VideoFrameRate/2)));
    fFeedStartxyz(i,3) = round(nanmean(mz(flightToFeederStart(i):flightToFeederStart(i)+trackData.VideoFrameRate/2)));
    
    fFeedEndxyz(i,1) = round(nanmean(mx(flightToFeederEnd(i):flightToFeederEnd(i)+trackData.VideoFrameRate/2)));
    fFeedEndxyz(i,2) = round(nanmean(my(flightToFeederEnd(i):flightToFeederEnd(i)+trackData.VideoFrameRate/2)));
    fFeedEndxyz(i,3) = round(nanmean(mz(flightToFeederEnd(i):flightToFeederEnd(i)+trackData.VideoFrameRate/2)));
    
    scatter3(fFeedStartxyz(i,1),fFeedStartxyz(i,2),fFeedStartxyz(i,3),100,'r','filled')
    hold on
    scatter3(fFeedEndxyz(i,1),fFeedEndxyz(i,2),fFeedEndxyz(i,3),100,'k','filled')
end
title(['Flights to feeder: ' batName ' ' dateSesh ' ' sessionType]);
xlabel('mm'); ylabel('mm'); zlabel('mm');
hold off

%plot flight paths of bat from feeder
flightFromFeederEnd = zeros(length(LT));
flightFromFeederStart = zeros(length(LT));
flightPathsFromFeeder = figure();
for i = 1:length(LT)
    flightFromFeederStart(i) = round(LT(i) * trackData.VideoFrameRate);
    flightFromFeederEnd(i) = flightFromFeederStart(i) + delay + (trackData.VideoFrameRate*preflightTime);
    %check that the end of the flight does not occur after session finishes
    if flightFromFeederEnd(i) > length(mx)
        flightFromFeederEnd(i) = length(mx)
    end
    plot3(mx(flightFromFeederStart(i):flightFromFeederEnd(i)),my(flightFromFeederStart(i):flightFromFeederEnd(i)),mz(flightFromFeederStart(i):flightFromFeederEnd(i)),'LineWidth',2,'color',CM(i,:))
    hold on
    %pause(0.1)
end
title(['Flights from feeder: ' batName ' ' dateSesh ' ' sessionType]);
xlabel('mm'); ylabel('mm'); zlabel('mm');
hold off

%k means cluster of flight trajectories into nclusters
%find pairs of start and endpoints with a high number of flights
rng(2) %control random number generation
kstart = kmeans(fFeedStartxyz,nclusters);
rng(2)
kend = kmeans(fFeedEndxyz,nclusters);
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

%plot clustered flights individually by each cluster
flightPathsClusterToFeederEach = figure();
jj = jet;
for traj = 1 : 10
    try
    if traj<ntrajectories + 1 % Only plot the first ntrajectories...
    subplot(round(ntrajectories/2),2,traj)
    for nf = allflights{ssf(traj)}
        plot3(mx(flightToFeederStart(nf):flightToFeederEnd(nf)),my(flightToFeederStart(nf):flightToFeederEnd(nf)),mz(flightToFeederStart(nf):flightToFeederEnd(nf)),'LineWidth',1,'Color',jj(traj*10,:))
        hold on
        scatter3(fFeedStartxyz(nf,1),fFeedStartxyz(nf,2),fFeedStartxyz(nf,3),100,'r','filled')
        hold on
        scatter3(fFeedEndxyz(nf,1),fFeedEndxyz(nf,2),fFeedEndxyz(nf,3),100,'k','filled')
        hold on
        %axis equal
        view(0,90)
        xlabel('mm'); ylabel('mm');
    end
    end
        clusterIndex{traj}(:) = allflights{ssf(traj)};

    catch
        disp(' no more flights...');
    end
end
sgtitle(['Flight Clusters to Feeder start(r)/stop(b): ' batName ' ' dateSesh ' ' sessionType]);
xlabel('mm'); ylabel('mm'); zlabel('mm');
hold off

%plot all the clusters in 1 figure in different colors
flightPathsClusterToFeederAll = figure();
jj = jet;
for traj = 1 : ntrajectories
    for nf = allflights{ssf(traj)}
        plot3(mx(flightToFeederStart(nf):flightToFeederEnd(nf)),my(flightToFeederStart(nf):flightToFeederEnd(nf)),mz(flightToFeederStart(nf):flightToFeederEnd(nf)),'LineWidth',1,'Color',jj(traj*10,:))
        hold on
        scatter3(fFeedStartxyz(nf,1),fFeedStartxyz(nf,2),fFeedStartxyz(nf,3),100,'r','filled')
        hold on
        scatter3(fFeedEndxyz(nf,1),fFeedEndxyz(nf,2),fFeedEndxyz(nf,3),100,'k','filled')
        hold on
    end
end
title(['Flight Clusters to Feeder start(r)/stop(b): ' batName ' ' dateSesh ' ' sessionType]);
%axis equal
view(0,90)
xlabel('mm'); ylabel('mm'); zlabel('mm');
hold off

flightFeedersStartStop.clusterIndex = clusterIndex;
flightFeedersStartStop.flightToFeederEnd = flightToFeederEnd;
flightFeedersStartStop.flightToFeederStart = flightToFeederStart;
flightFeedersStartStop.flightFromFeederEnd = flightFromFeederEnd;
flightFeedersStartStop.flightFromFeederStart = flightFromFeederStart;
flightFeedersStartStop.fFeedStartxyz = fFeedStartxyz;
flightFeedersStartStop.fFeedEndxyz = fFeedEndxyz;
catch
    
    flightFeedersStartStop.clusterIndex = [];
    flightPathsFromFeeder = figure();
    flightPathsClusterToFeederEach = figure();
    flightPathsClusterToFeederAll = figure();
    flightFeedersStartStop.flightToFeederEnd = [];
    flightFeedersStartStop.flightToFeederStart = [];
    flightFeedersStartStop.flightFromFeederEnd = [];
    flightFeedersStartStop.flightFromFeederStart = [];
    flightFeedersStartStop.fFeedStartxyz = [];
    flightFeedersStartStop.fFeedEndxyz = [];
    
end


