function [flightPathsAll,flightPathsStartStop, flightPaths, flightPathsClusterEach, flightPathsClusterAll] = ImBat_plotFlights(trackData)

global batName dateSesh sessionType

nclusters = 4; %number of clusters for kmeans clustering of flight trajectories
ntrajectories = 6; %number of output trajectories from kmeans that you want to look at 

%find ttl pulses for synching
%event_ttls = trackData.AnalogSignals(:,2); %from motion data
%[R,LT,UT,LL,UL] = risetime(event_ttls,trackData.VideoFrameRate); %find times of ttl pulses in SECONDS

%might want to use raw c3d to keep track of each marker (this could be
%useful for when not all 3 markers are present, which is often the case
%when the animal is at a wall)
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

%threshold based on speed
Vx = gradient(mx, 1/trackData.VideoFrameRate);
Vy = gradient(my, 1/trackData.VideoFrameRate);
Vz = gradient(mz, 1/trackData.VideoFrameRate);
batSpeed = sqrt(Vx.^2 + Vy.^2 + Vz.^2)/1000; %in m/s

nonflying = find(batSpeed < 1.5);
toofast = find(batSpeed > 10);

%set low speed points to zero
mx([nonflying; toofast]) = nan;
my([nonflying; toofast]) = nan;
mz([nonflying; toofast]) = nan;

trajectories_continuous(1,:) = mx';
trajectories_continuous(2,:) = my';
trajectories_continuous(3,:) = mz';

%plot all flights in black
flightPathsAll = figure();
plot3(mx,my,mz,'LineWidth',2,'Color','k')
hold on
% % modify labels for tick marks
% xticks = get(gca,'xtick');
% yticks = get(gca,'ytick');
% zticks = get(gca,'ztick');
% scaling  = 0.1; %1.1um per pixel
% newlabelsX = arrayfun(@(ax) sprintf('%g', scaling * ax), xticks, 'un', 0);
% newlabelsY = arrayfun(@(ay) sprintf('%g', scaling * ay), yticks, 'un', 0);
% newlabelsZ = arrayfun(@(az) sprintf('%g', scaling * az), zticks, 'un', 0);
% set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY,'zticklabel',newlabelsZ);
title(['All flights: ' batName ' ' dateSesh ' ' sessionType]);
xlabel('mm'); ylabel('mm'); zlabel('mm');
view(0,90)
hold off

%splice out individual flights
batspeed = batSpeed;
batspeed(nonflying) = nan;
batspeed(toofast) = nan;

bflying=~isnan(batspeed);

%for each data point, sum up the next 1s of data
allsums = [];
for bf = 1 : size(bflying,1)-trackData.VideoFrameRate
    allsums(bf) = sum(bflying(bf:bf+trackData.VideoFrameRate/2));
end

[R,rLT,rUT,rLL,rUL] = risetime(allsums);
[F,fLT,fUT,fLL,fUL] = falltime(allsums);
flight_starts = round(rLT);
flight_ends = round(fLT);

%cut out flights and save
flightPathsStartStop = figure();
if size(R,2) > 0
    CM = jet(size(R,2));
    for nf = 1 : size(R,2)
        hold on
        plot3(mx(flight_starts(nf):flight_ends(nf)),my(flight_starts(nf):flight_ends(nf)),mz(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',CM(nf,:))
        hold on
        
        fstartxyz(nf,1) = round(nanmean(mx(flight_starts(nf):flight_starts(nf)+trackData.VideoFrameRate/2)));
        fstartxyz(nf,2) = round(nanmean(my(flight_starts(nf):flight_starts(nf)+trackData.VideoFrameRate/2)));
        fstartxyz(nf,3) = round(nanmean(mz(flight_starts(nf):flight_starts(nf)+trackData.VideoFrameRate/2)));
        
        fendxyz(nf,1) = round(nanmean(mx(flight_ends(nf):flight_ends(nf)+trackData.VideoFrameRate/2)));
        fendxyz(nf,2) = round(nanmean(my(flight_ends(nf):flight_ends(nf)+trackData.VideoFrameRate/2)));
        fendxyz(nf,3) = round(nanmean(mz(flight_ends(nf):flight_ends(nf)+trackData.VideoFrameRate/2)));
        
        scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),100,'r','filled')
        hold on
        scatter3(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),100,'k','filled')
        %pause
    end
else
    fstartxyz(1,1) = (0);
    fstartxyz(1,2) = (0);
    fstartxyz(1,3) = (0);
    
    fendxyz(1,1) = (0);
    fendxyz(1,2) = (0);
    fendxyz(1,3) = (0);    
end
title(['All flights start(r)/stop(b): ' batName ' ' dateSesh ' ' sessionType]);
xlabel('mm'); ylabel('mm'); zlabel('mm');
hold off

%k means cluster of flight trajectories into nclusters
%find pairs of start and endpoints with a high number of flights
rng(2) %control random number generation
kstart = kmeans(fstartxyz,nclusters);
rng(2)
kend = kmeans(fendxyz,nclusters);
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
%plot clustered flights
flightPathsClusterEach = figure();
jj = jet;
for traj = 1 : 10
    try
    if traj<ntrajectories + 1 % Only plot the first ntrajectories...
    subplot(3,2,traj)
    for nf = allflights{ssf(traj)}
        plot3(mx(flight_starts(nf):flight_ends(nf)),my(flight_starts(nf):flight_ends(nf)),mz(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',jj(traj*10,:))
        hold on
        scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),100,'r','filled')
        hold on
        scatter3(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),100,'k','filled')
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
sgtitle(['Flight Clusters start(r)/stop(b): ' batName ' ' dateSesh ' ' sessionType]);
xlabel('mm'); ylabel('mm'); zlabel('mm');
hold off

%plot all the clusters in 1 figure in different colors
flightPathsClusterAll = figure();
jj = jet;
for traj = 1 : ntrajectories
    for nf = allflights{ssf(traj)}
        plot3(mx(flight_starts(nf):flight_ends(nf)),my(flight_starts(nf):flight_ends(nf)),mz(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',jj(traj*10,:))
        hold on
        scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),100,'r','filled')
        hold on
        scatter3(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),100,'k','filled')
        hold on
    end
end
title(['Flight Clusters start(r)/stop(b): ' batName ' ' dateSesh ' ' sessionType]);
%axis equal
view(0,90)
xlabel('mm'); ylabel('mm'); zlabel('mm');
hold off


flightPaths.flight_starts_idx = flight_starts;
flightPaths.flight_ends_idx = flight_ends;
flightPaths.flight_ends_xyz = fendxyz;
flightPaths.flight_starts_xyz = fstartxyz;
flightPaths.trajectories_continuous = trajectories_continuous; 
flightPaths.batSpeed = batSpeed; 
flightPaths.clusterIndex = clusterIndex;

