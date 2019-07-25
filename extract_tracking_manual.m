idx = Markers == 0;
Markers(idx) = NaN;
avgMarker = squeeze(nanmean(Markers,2));
figure;
scatter3(avgMarker(:,1),avgMarker(:,2),avgMarker(:,3))

%%
figure
for i = 1:15
scatter3(Markers(:,i,1),Markers(:,i,2),Markers(:,i,3))
hold on
pause 
end

%%
event_ttls = AnalogSignals(:,1); %trial data
[R,LT,UT,LL,UL] = risetime(event_ttls,VideoFrameRate);
markerSet = Markers;
figure

for i = 1:length(LT)
    flightEnd = round(LT(i) * 120);
    flightStart = flightEnd - 240;
    scatter3(markerSet(flightStart:flightEnd,1,1),markerSet(flightStart:flightEnd,1,2),markerSet(flightStart:flightEnd,1,3))
    hold on
    pause(0.1)
end
