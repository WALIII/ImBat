trackData = load('tracking data file.mat')

rewardSignal = trackData.AnalogSignals(:,1);
[rewardPks,rewardLocs] = findpeaks(rewardSignal,'MinPeakProminence',2);


%plot ttl times of reward landing
figure(); 
plot(rewardSignal);
hold on;
text(rewardLocs+.02,rewardPks,num2str((1:numel(rewardPks))'));
