function ImBat_FindEvents(neuron, out)


% sort data to correct for weird timstamps
[~, out.video_times_true] =  sort(out.video_times);
%Cx = neuron.S(:,out.video_times);

% Movement events:
figure();
D = diff(nanmean(out.Location2(:,1:3),2));
plot(out.Location_time(1:end-1),(abs(D)))


Cx = zscore(smooth(mean(full(neuron.S(1:50,out.video_times(1:end-1))),1),200));
[CxA CxB] = findpeaks(Cx,out.video_times_true(1:end-1),'MinPeakProminence',4);

% Neural events:
figure();
hold on;
plot(out.video_times(1:end-1),Cx);
plot(CxB,CxA);
ylim([0 5]);



% Reactivatons 

