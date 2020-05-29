close all;
batName = 'Gio';
dateSesh = '200521';
sessionType = 'OnOffOn';

duration = length(Markers(:,1,1));

mx = Markers(:,1,1);%trackData.Markers(1:length(trackData.Markers(:,1,1))/2,1,1);%trackData.Markers(length(trackData.Markers(:,1,1))/2+1:length(trackData.Markers(:,1,1)),1,1);
my = Markers(:,1,2);%trackData.Markers(1:length(trackData.Markers(:,1,1))/2,1,2);%trackData.Markers(length(trackData.Markers(:,1,1))/2+1:length(trackData.Markers(:,1,1)),1,2);
mz = Markers(:,1,3);%trackData.Markers(1:length(trackData.Markers(:,1,1))/2,1,3);%trackData.Markers(length(trackData.Markers(:,1,1))/2+1:length(trackData.Markers(:,1,1)),1,3);

%set zeros to nan
mx(find(mx == 0)) = nan;
my(find(my == 0)) = nan;
mz(find(mz == 0)) = nan;

mxFull = fillmissing(mx,'nearest');
myFull = fillmissing(my,'nearest');
mzFull = fillmissing(mz,'nearest');

rgbOnOffOn = figure();
sgtitle([batName ' ' dateSesh ' ' sessionType]);
subplot(2,2,1)
plot3(mxFull(1:(duration/3),1,1),myFull(1:(duration/3),1),mzFull(1:(duration/3),1),'r')
hold on;
plot3(mxFull((duration/3):(duration*2/3),1,1),myFull((duration/3):(duration*2/3),1),mzFull((duration/3):(duration*2/3),1),'g')
hold on
plot3(mxFull((duration*2/3):end,1),myFull((duration*2/3):end,1),mzFull((duration*2/3):end,1),'b')
view(0,90);
xlim([-3000 3000]);
ylim([-3000 3000]);
title('60 min session');

subplot(2,2,2)
plot3(mxFull(1:(duration/3),1,1),myFull(1:(duration/3),1),mzFull(1:(duration/3),1),'r')
view(0,90);
xlim([-3000 3000]);
ylim([-3000 3000]);
title('First 20 min');
subplot(2,2,3)
plot3(mxFull((duration/3):(duration*2/3),1),myFull((duration/3):(duration*2/3),1),mzFull((duration/3):(duration*2/3),1),'g')
view(0,90);
xlim([-3000 3000]);
ylim([-3000 3000]);
title('Middle 20 min');
subplot(2,2,4)
plot3(mxFull((duration*2/3):end,1),myFull((duration*2/3):end,1),mzFull((duration*2/3):end,1),'b')
view(0,90);
xlim([-3000 3000]);
ylim([-3000 3000]);
title('Last 20 min');


saveas(rgbOnOffOn,['/Users/periscope/Desktop/RGB_flights_onOffOn/' batName '_' dateSesh '_' sessionType '.tif']);