
fs = 120;
batName = 'Gen';
dateSesh = '200530';
sessionType = 'fly-1';
[rewardR,rewardLT,rewardUT] = risetime(out.RewardVector);
rewardYes_idx = find(bhv_data.trials(:,9)==1);
rewardNo_idx = find(bhv_data.trials(:,9)==0);
rewardYes_times = rewardLT(rewardYes_idx);
rewardNo_times = rewardLT(rewardNo_idx);

figure();
subplot(1,2,1)
for t = 1:length(rewardYes_times)
[minEndYes(t),minEndYesIdx(t)] = min(abs(rewardYes_times(t)-flightPaths.flight_ends_idx));
plot3(flightPaths.trajectoriesContinous(1,flightPaths.flight_starts_idx(minEndYesIdx(t)):flightPaths.flight_ends_idx(minEndYesIdx(t))),flightPaths.trajectoriesContinous(2,flightPaths.flight_starts_idx(minEndYesIdx(t)):flightPaths.flight_ends_idx(minEndYesIdx(t))),flightPaths.trajectoriesContinous(3,flightPaths.flight_starts_idx(minEndYesIdx(t)):flightPaths.flight_ends_idx(minEndYesIdx(t))));
hold on;
%scatter3(flightPaths.flight_starts_xyz(1,1,t),flightPaths.flight_starts_xyz(2,1,t),flightPaths.flight_starts_xyz(3,1,t),50,'r','filled')
%scatter3(flightPaths.flight_ends_xyz(1,1,t),flightPaths.flight_ends_xyz(2,1,t),flightPaths.flight_ends_xyz(3,1,t),50,'r','filled')
end
xlim([-3 3])
ylim([-3 3])
view(0,90)
xlabel('m'); ylabel('m');
title(['Rewarded flights: ' batName ' ' dateSesh ' ' sessionType]);

subplot(1,2,2)
for t = 1:length(rewardNo_times)
[minEndNo(t),minEndNoIdx(t)] = min(abs(rewardNo_times(t)-flightPaths.flight_ends_idx));
plot3(flightPaths.trajectoriesContinous(1,flightPaths.flight_starts_idx(minEndNoIdx(t)):flightPaths.flight_ends_idx(minEndNoIdx(t))),flightPaths.trajectoriesContinous(2,flightPaths.flight_starts_idx(minEndNoIdx(t)):flightPaths.flight_ends_idx(minEndNoIdx(t))),flightPaths.trajectoriesContinous(3,flightPaths.flight_starts_idx(minEndNoIdx(t)):flightPaths.flight_ends_idx(minEndNoIdx(t))));
hold on;
end
xlim([-3 3])
ylim([-3 3])
view(0,90)
xlabel('m'); ylabel('m');
title(['Catch trials: ' batName ' ' dateSesh ' ' sessionType]);

for i = 1:length(flightPaths.clusterIndex{2})
targetPos(i,:) = flightPaths.flight_ends_xyz(flightPaths.clusterIndex{2}(i),:);
targetTime(i) = flightPaths.flight_ends_idx(flightPaths.clusterIndex{2}(i));
[minStart(i),minStartIdx(i)] = min(abs(rewardYes_times-targetTime(i)));
%flightPaths.pos(:,flightPaths.flight_ends_idx(i),i)    
end
feeder2 = mean(targetPos);

