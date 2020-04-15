function [flightPaths] = ImBat_plotFlights(trackData,varargin)


nclusters = 6; %nIumber of clusters for kmeans clustering of flight trajectories
ntrajectories = 6; %number of output trajectories from kmeans that you want to look at


batName = [];
dateSesh = [];
sessionType = [];
saveFlag = 0;
loadFlag = 0;

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
        case 'loadflag'
            loadFlag = varargin{i+1};
        case 'analysisfolder'
            analysis_Folder = varargin{i+1};
    end
end

%labels for loading and saving data if running independent fromImBat_analyze
if loadFlag == 1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
    folderName = extractBefore(pwd,batName);
    %trackData = load([folderName label '_track.mat']);
end
%find ttl pulses for synching
%event_ttls = trackData.AnalogSignals(:,2); %from motion data
%[R,LT,UT,LL,UL] = risetime(event_ttls,trackData.VideoFrameRate); %find times of ttl pulses in SECONDS

%might want to use raw c3d to keep track of each marker (this could be
%useful for when not all 3 markers are present, which is often the case
%when the animal is at a wall)
idx = trackData.Markers == 0;
trackData.Markers(idx) = NaN;
avgMarker = squeeze(nanmean(trackData.Markers,2));
mx=avgMarker(:,1);
my=avgMarker(:,2);
mz=avgMarker(:,3);
%if using bat cluster
mx = trackData.Markers(:,1,1);%trackData.Markers(1:length(trackData.Markers(:,1,1))/2,1,1);%trackData.Markers(length(trackData.Markers(:,1,1))/2+1:length(trackData.Markers(:,1,1)),1,1);
my = trackData.Markers(:,1,2);%trackData.Markers(1:length(trackData.Markers(:,1,1))/2,1,2);%trackData.Markers(length(trackData.Markers(:,1,1))/2+1:length(trackData.Markers(:,1,1)),1,2);
mz = trackData.Markers(:,1,3);%trackData.Markers(1:length(trackData.Markers(:,1,1))/2,1,3);%trackData.Markers(length(trackData.Markers(:,1,1))/2+1:length(trackData.Markers(:,1,1)),1,3);

%set zeros to nan
mx(find(mx == 0)) = nan;
my(find(my == 0)) = nan;
mz(find(mz == 0)) = nan;

trajPartial(1,:) = mx';
trajPartial(2,:) = my';
trajPartial(3,:) = mz';

%mxNan = find(mx == nan);
mxNan = find(isnan(mx));

mxFull = fillmissing(mx,'nearest');
myFull = fillmissing(my,'nearest');
mzFull = fillmissing(mz,'nearest');

%threshold based on speed
Vx = gradient(mxFull, 1/trackData.VideoFrameRate);
Vy = gradient(myFull, 1/trackData.VideoFrameRate);
Vz = gradient(mzFull, 1/trackData.VideoFrameRate);
batSpeed = sqrt(Vx.^2 + Vy.^2 + Vz.^2)/1000; %in m/s

nonflying = find(batSpeed < 1.5);
toofast = find(batSpeed > 35);

%set low speed points to zero
mxFull([nonflying]) = nan;%mx(nonflying) = nan; %mx([nonflying; toofast]) = nan;
myFull([nonflying]) = nan;%my(nonflying) = nan;%my([nonflying; toofast]) = nan;
mzFull([nonflying]) = nan;%mz(nonflying) = nan;%mz([nonflying; toofast]) = nan;



trajectories_continuous(1,:) = mxFull';
trajectories_continuous(2,:) = myFull';
trajectories_continuous(3,:) = mzFull';

%% plot all flights in black
plotFlightPathsAll = figure();
plot3(mx,my,mz,'LineWidth',2,'Color','k')
%scatter3(mxFull,myFull,mzFull,'.','k')
hold on

% col  = jet(100);
% for i = 2:length(trackData.Markers(:,1,:))
% plot3(mxFull(i), myFull(i), mzFull(i), 'color',col(vecnorm([mxFull(i),myFull(i),mzFull(i) - mxFull(i-1),myFull(i-1),mzFull(i-1)], 2, 2)))
% hold on
% end

% % modify labels for tick marks
view(0,90)
xlim([-3000 3000])
ylim([-3000 3000])
xticks = get(gca,'xtick');
yticks = get(gca,'ytick');
newlabelsX = arrayfun(@(ax) sprintf('%g', ax/1000), xticks, 'un', 0);
newlabelsY = arrayfun(@(ay) sprintf('%g', ay/1000), yticks, 'un', 0);
set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
title(['All flights: ' batName ' ' dateSesh ' ' sessionType]);
xlabel('m'); ylabel('m');
hold off

%% splice out individual flights
batspeed = batSpeed;
batspeed(nonflying) = nan;
%batspeed(toofast) = nan;

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
if length(flight_starts) ~= length(flight_ends)
    flight_ends = [flight_ends size(allsums,2)];
end

%% cut out flights and save
plotFlightPathsStartStop = figure();
if size(R,2) > 0
    CM = jet(size(R,2));
    for nf = 1 : size(R,2)
        hold on
        plot3(mxFull(flight_starts(nf):flight_ends(nf)),myFull(flight_starts(nf):flight_ends(nf)),mzFull(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',CM(nf,:))
        hold on
        
        fstartxyz(nf,1) = round(nanmean(mxFull(flight_starts(nf):flight_starts(nf)+trackData.VideoFrameRate/2)));
        fstartxyz(nf,2) = round(nanmean(myFull(flight_starts(nf):flight_starts(nf)+trackData.VideoFrameRate/2)));
        fstartxyz(nf,3) = round(nanmean(mzFull(flight_starts(nf):flight_starts(nf)+trackData.VideoFrameRate/2)));
        
        fendxyz(nf,1) = round(nanmean(mxFull(flight_ends(nf):flight_ends(nf)+trackData.VideoFrameRate/2)));
        fendxyz(nf,2) = round(nanmean(myFull(flight_ends(nf):flight_ends(nf)+trackData.VideoFrameRate/2)));
        fendxyz(nf,3) = round(nanmean(mzFull(flight_ends(nf):flight_ends(nf)+trackData.VideoFrameRate/2)));
        
        scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),50,'r','filled')
        hold on
        scatter3(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),50,'k','filled')
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
% modify labels for tick marks
view(0,90)
xlim([-3000 3000])
ylim([-3000 3000])
xticks = get(gca,'xtick');
yticks = get(gca,'ytick');
newlabelsX = arrayfun(@(ax) sprintf('%g', ax/1000), xticks, 'un', 0);
newlabelsY = arrayfun(@(ay) sprintf('%g', ay/1000), yticks, 'un', 0);
set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
xlabel('m'); ylabel('m');

hold off
%try
%% k means cluster of flight trajectories into nclusters
%find pairs of start and endpoints with a high number of flights
rng(2) %control random number generation
try
    kstart = kmeans(fstartxyz,nclusters);
catch
    kstart = kmeans(fstartxyz,round(nclusters/2));
end
rng(2)
try
    kend = kmeans(fendxyz,nclusters);
catch
    kend = kmeans(fendxyz,round(nclusters/2));
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
%% plot clustered flights
plotFlightPathsClusterEach = figure();
sgtitle(['Flight Clusters start(r)/stop(b): ' batName ' ' dateSesh ' ' sessionType]);
jj = jet(ntrajectories);
for traj = 1 : 10
    try
        if traj<ntrajectories + 1 % Only plot the first ntrajectories...
            subplot(3,2,traj)
            for nf = allflights{ssf(traj)}
                plot(mxFull(flight_starts(nf):flight_ends(nf)),myFull(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',jj(traj,:))
                %plot3(mxFull(flight_starts(nf):flight_ends(nf)),myFull(flight_starts(nf):flight_ends(nf)),mzFull(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',jj(traj,:))
                hold on
                scatter(fstartxyz(nf,1),fstartxyz(nf,2),50,'r','filled')
                %scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),50,'r','filled')
                hold on
                scatter(fendxyz(nf,1),fendxyz(nf,2),50,'k','filled')
                %scatter3(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),50,'k','filled')
                hold on
                %axis equal
                % % modify labels for tick marks
                %view(0,90)
                xlim([-3000 3000])
                ylim([-3000 3000])
                xticks = get(gca,'xtick');
                yticks = get(gca,'ytick');
                newlabelsX = arrayfun(@(ax) sprintf('%g', ax/1000), xticks, 'un', 0);
                newlabelsY = arrayfun(@(ay) sprintf('%g', ay/1000), yticks, 'un', 0);
                set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
                xlabel('m'); ylabel('m');
                title(['Cluster ' num2str(traj) ':' num2str(length(allflights{ssf(traj)})) ' flights']);
            end
        end
        clusterIndex{traj}(:) = allflights{ssf(traj)};
        
    catch
        disp(' no more flights...');
    end
end
hold off

%% plot all the clusters in 1 figure in different colors
plotFlightPathsClusterAll = figure();
jj = jet(ntrajectories);
for traj = 1 : ntrajectories
    for nf = allflights{ssf(traj)}
        %plot(mxFull(flight_starts(nf):flight_ends(nf)),myFull(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',jj(traj,:))
        plot3(mxFull(flight_starts(nf):flight_ends(nf)),myFull(flight_starts(nf):flight_ends(nf)),mzFull(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',jj(traj,:))
        hold on
        scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),50,'r','filled')
        hold on
        %scatter(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),50,'k','filled')
        scatter3(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),50,'k','filled')
        hold on
    end
end
title(['Flight Clusters start(r)/stop(b): ' batName ' ' dateSesh ' ' sessionType]);
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

        hold off%%   
flightPaths.clusterIndex = clusterIndex;
% catch
% flightPathsClusterEach = figure();
% flightPathsClusterAll = figure();
% flightPaths.clusterIndex = [];
% end
flightPaths.flight_starts_idx = flight_starts;
flightPaths.flight_ends_idx = flight_ends;
flightPaths.flight_ends_xyz = fendxyz;
flightPaths.flight_starts_xyz = fstartxyz;
flightPaths.trajectories_continuous = trajectories_continuous;
flightPaths.batSpeed = batSpeed;
flightPaths.flightPathsAll = plotFlightPathsAll;
flightPaths.flightPathsClusterEach = plotFlightPathsClusterEach;
flightPaths.flightPathsClusterAll = plotFlightPathsClusterAll;
flightPaths.flightPathsStartStop = plotFlightPathsStartStop;


if saveFlag == 1
    mkdir([analysis_Folder '\flights']);
    %set(findall(plotFlightPathsAll,'-property','FontSize'),'FontSize',20);
    saveas(plotFlightPathsAll,[pwd '\' analysis_Folder '\flights\' batName '_' dateSesh '_' sessionType '_flightPathsAll_full.tif']);
    savefig(plotFlightPathsAll,[pwd '\' analysis_Folder '\flights\' batName '_' dateSesh '_' sessionType '_flightPathsAll_full.fig']);
    saveas(plotFlightPathsAll,[pwd '\' analysis_Folder '\flights\' batName '_' dateSesh '_' sessionType '_flightPathsAll_full.svg']);
    %set(findall(plotFlightPathsClusterEach,'-property','FontSize'),'FontSize',20);
    saveas(plotFlightPathsClusterEach,[pwd '\' analysis_Folder '\flights\' batName '_' dateSesh '_' sessionType '_flightPathsClusterEach_full.tif']);
    savefig(plotFlightPathsClusterEach,[pwd '\' analysis_Folder '\flights\' batName '_' dateSesh '_' sessionType '_flightPathsClusterEach_full.fig']);
    saveas(plotFlightPathsClusterEach,[pwd '\' analysis_Folder '\flights\' batName '_' dateSesh '_' sessionType '_flightPathsClusterEach_full.svg']);
    %set(findall(plotFlightPathsClusterAll,'-property','FontSize'),'FontSize',20);
    saveas(plotFlightPathsClusterAll,[pwd '\' analysis_Folder '\flights\' batName '_' dateSesh '_' sessionType '_flightPathsClusterAll_full.tif']);
    savefig(plotFlightPathsClusterAll,[pwd '\' analysis_Folder '\flights\' batName '_' dateSesh '_' sessionType '_flightPathsClusterAll_full.fig']);
    saveas(plotFlightPathsClusterAll,[pwd '\' analysis_Folder '\flights\' batName '_' dateSesh '_' sessionType '_flightPathsClusterAll_full.svg']);
    %set(findall(plotFlightPathsStartStop,'-property','FontSize'),'FontSize',20);
    saveas(plotFlightPathsStartStop,[pwd '\' analysis_Folder '\flights\' batName '_' dateSesh '_' sessionType '_flightPathsStartStop_full.tif']);
    savefig(plotFlightPathsStartStop,[pwd '\' analysis_Folder '\flights\' batName '_' dateSesh '_' sessionType '_flightPathsStartStop_full.fig']);
    saveas(plotFlightPathsStartStop,[pwd '\' analysis_Folder '\flights\' batName '_' dateSesh '_' sessionType '_flightPathsStartStop_full.svg']);

    save([pwd '\' analysis_Folder '\' batName '_' dateSesh '_' sessionType '_flightPaths.mat'],'flightPaths')
end

