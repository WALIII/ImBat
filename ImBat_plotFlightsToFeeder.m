function [flightPathsToFeeder, flightPathsFromFeeder, flightToFeederEnd, flightToFeederStart, flightFromFeederEnd, flightFromFeederStart]...
    = ImBat_plotFlightsToFeeder(trackData)

global batName dateSesh sessionType

%length of time for flight up until reaching feeder (s)
preflightTime = 2;
%length of time to wait before bat leaves feeder for FlightFromFeeder track %(s)
delay = 3;

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
for i = 1:length(LT)
    flightToFeederEnd(i) = round(LT(i) * trackData.VideoFrameRate);
    flightToFeederStart(i) = flightToFeederEnd(i) - (trackData.VideoFrameRate*preflightTime); %minus the preflight # of seconds
    %scatter3(markerSet(flightStart:flightEnd,1,1),markerSet(flightStart:flightEnd,1,2),markerSet(flightStart:flightEnd,1,3))
    plot3(mx(flightToFeederStart(i):flightToFeederEnd(i)),my(flightToFeederStart(i):flightToFeederEnd(i)),mz(flightToFeederStart(i):flightToFeederEnd(i)),'LineWidth',2,'color',CM(i,:))
    hold on
    %pause(0.1)
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
    plot3(mx(flightFromFeederStart(i):flightFromFeederEnd(i)),my(flightFromFeederStart(i):flightFromFeederEnd(i)),mz(flightFromFeederStart(i):flightFromFeederEnd(i)),'LineWidth',2,'color',CM(i,:))
    hold on
    %pause(0.1)
end
title(['Flights from feeder: ' batName ' ' dateSesh ' ' sessionType]);
xlabel('mm'); ylabel('mm'); zlabel('mm');
hold off