function [flightPathsAll,flightPathsStartStop, trajectories_continuous, flight_starts, flight_ends, fstartxyz, fendxyz, batSpeed] = ImBat_plotFlights(trackData)

global batName dateSesh sessionType

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
        
        scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),300,'g','filled')
        hold on
        scatter3(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),300,'b','filled')
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
title(['All flights start(g)/stop(b): ' batName ' ' dateSesh ' ' sessionType]);
xlabel('mm'); ylabel('mm'); zlabel('mm');
hold off

