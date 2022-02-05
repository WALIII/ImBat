function [out] = ImBat_Goal_Distance(flightPaths,FlightAlignedROI,reward_location,Reward_Number);

% WAL3
% Nov 26, 2021

if isempty(Reward_Number)
Reward_Number =1;
end
% Preprocessin:
disp('aligning data');
[S2save,F2save,other] =  ImBat_Spikes(FlightAlignedROI);


% Calculate the distance and angle to the goal 


% User params:
numROIs = size(S2save,2);
distance_smooth = 10;
angle_smooth = 10;
spike_thresh = 0.05;
reward_toggle = 0; % look at second reward location if '1'
Euclid_distance = 0; % use one if you want euclid distance to goal,else distance along flight path. 
velocity_thresh = 0;%1.5;
distance_thresh = 400;
col = {'r','g'};


if isempty(reward_location)
    
% look at reward locations first..
figure();
hold on;
 plot(flightPaths.AllFlights(:,1),flightPaths.AllFlights(:,2),'-')
 plot(flightPaths.AllFlights(flightPaths.RewardIdx,1),flightPaths.AllFlights(flightPaths.RewardIdx,2),'*')

 clear Reward_Locations
 Reward_Locations(:,1) = flightPaths.AllFlights(flightPaths.RewardIdx,1);
 Reward_Locations(:,2) = flightPaths.AllFlights(flightPaths.RewardIdx,2);

 [idx,C] = kmeans(Reward_Locations,2); % cluster reward triggers 
 
 figure();
hold on;
 plot(flightPaths.AllFlights(:,1),flightPaths.AllFlights(:,2),'-');
 for i = 1:2
 plot(flightPaths.AllFlights(flightPaths.RewardIdx(idx==i),1),flightPaths.AllFlights(flightPaths.RewardIdx(idx==i),2),'*','color',col{i} )
 
 M(1) = median(flightPaths.AllFlights(flightPaths.RewardIdx(idx==i),1));
 M(2) = median(flightPaths.AllFlights(flightPaths.RewardIdx(idx==i),2));
  plot(M(1),M(2),'o','MarkerSize',50,'color',col{i})
  reward_loc{i}.M = M;
 end
 title( 'Clustered Reward Locations');
 
% Check if the top reward is designated as 'reward location 1';
if Reward_Number == 1;
if reward_loc{1}.M(2)>reward_loc{2}.M(2)
else
   reward_loc_temp{1} = reward_loc{2};
   reward_loc_temp{2} = reward_loc{1};
   reward_loc = reward_loc_temp;
end
else
    % flip the reward if you have the 'Reward_Number' ==2
  if reward_loc{1}.M(2)<reward_loc{2}.M(2)
else
   reward_loc_temp{2} = reward_loc{1};
   reward_loc_temp{1} = reward_loc{2};
   reward_loc = reward_loc_temp;
  end
end
   
else
     reward_loc = reward_location;
end
out.reward_location = reward_loc;

% remove indexes where velocity is too small (will mess up angle measurements at landing sites)).
ind2rm = find(other.T_spd<velocity_thresh);
F2save(ind2rm,:) = [];
S2save(ind2rm,:) = [];
other.T_dist(ind2rm) = [];
ind2rm2 = find(other.T_dist<distance_thresh);
F2save(ind2rm2,:) = [];
S2save(ind2rm2,:) = [];
other.T_dist(ind2rm2) = [];
%%% Velocity
% clear x1 x2 y1 y2 z1 z2 s spd
x1 = F2save(2:end,1);
y1 = F2save(2:end,2);
z1 = F2save(2:end,3);
% 
% S=sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);
% S(S>40) = 0;
spd = other.T_spd;

theta = zeros(1, size(x1,1));


%%% angle (azmuthal) offset by 60 samples.
for i = 1:size(x1,1)-60
    theta(i) = atan2(y1(i+60)-y1(i),x1(i+60)-x1(i));
end
theta = theta';


% Goal calculations
% compute angle to the goal

for i = 1:size(x1,1)-60
   theta_reward1(i) = atan2(reward_loc{1}.M(2)-y1(i),reward_loc{1}.M(1)-x1(i));
    theta_reward2(i) = atan2(reward_loc{2}.M(2)-y1(i),reward_loc{2}.M(1)-x1(i));
%      theta_reward1(i) = atan2(y1(i)-reward_loc{1}.M(2),x1(i)-reward_loc{1}.M(1));
%    theta_reward2(i) = atan2(y1(i)-reward_loc{2}.M(2),x1(i)-reward_loc{2}.M(1));
end



% calculte reward distance:
if Euclid_distance == 1;
distance_reward1 = sqrt( (reward_loc{1}.M(2)-y1).^2 + (reward_loc{1}.M(1)-x1).^2 );
distance_reward2 = sqrt( (reward_loc{2}.M(2)-y1).^2 + (reward_loc{1}.M(2)-x1).^2 );
else
    disp('warning: using end of flight instead of euclid distance to goal');
    distance_reward1 = other.T_dist;
    distance_reward2 = other.T_dist;
end

% units of degrees
h1 =  rad2deg(theta)';
D1 = rad2deg(theta_reward1);  % location 1
D2 = rad2deg(theta_reward2); % location 2


% subtract current heading 
D1 = D1-h1(1:size(D1,2)); % normalize to headding 
D2 = D2-h1(1:size(D1,2));

% toggle
if reward_toggle ==1;
D1 = D2;
end

% offsets
D1(D1<-180) = D1(D1<-180)+360;
D1(D1>180) = D1(D1>180)-360;
D2(D2<-180) = D2(D2<-180)+360;
D2(D2>180) = D2(D2>180)-360;

% plot the distance vs theta distributions..
figure(); plot(other.T_dist(1:size(theta_reward1,2)), D1,'.');

% Bin spikes..
figure();
for i = 1:numROIs
ROI = i;

ind2try = find(S2save(:,ROI)>spike_thresh);
ind2try(ind2try>size(D1,2)) = []; % get rid of edge..

% Calculate binned spikes vs angle or distance
angle_spikes = D1(ind2try);
distance_spikes = distance_reward1(ind2try);

% ANGLE CALCULTAION
A1 = histcounts(angle_spikes,'Normalization','probability','Binedges',[-180:30:180]);
B1 = histcounts(D1+0.1,'Normalization','probability','Binedges',[-180:30:180]);
norm_angle = A1./B1;
norm_angle(isnan(norm_angle)) = 0;
bar(norm_angle)
xticks([1 18 36])
ylabel('ratio relative to chance')
xlabel('Goal angle');
xticklabels({'-180','0','180'})
title('unclustered, sorted by clustered');
n(i,:) = mat2gray(smooth(norm_angle,angle_smooth)); 

comp_vect = smooth(norm_angle,angle_smooth);
rl_pk_angle(i) = max(comp_vect)-min(comp_vect); % angle real peak ( to compare to shuff)
clear A1
% pause();
 
 % DISTANCE CALCULATION
[dA1 edgesdA1] = histcounts(distance_spikes,'Normalization','probability','BinEdges',[distance_thresh:100:10000]);
[dB1 edgesdB1] = histcounts(distance_reward1,'Normalization','probability','BinEdges',[distance_thresh:100:10000]);
norm_distance = dA1./dB1;
norm_distance(isnan(norm_distance)) = 0;

bar(norm_distance)

n2(i,:) = mat2gray(smooth(norm_distance,distance_smooth)); 
comp_vect2 = smooth(norm_distance,distance_smooth);
rl_pk_distance(i) = max(comp_vect2)-min(comp_vect2); % distance real peak

% xticks([1 18 36])
ylabel('ratio relative to chance')
xlabel('Goal angle');
xticklabels({'-180','0','180'})
title('unclustered, sorted by clustered');
%  pause();
end


[~, g1] = max(n');
 [~,st1] = sort(g1);
figure(); imagesc((n(st1,:)));
title('angle snake');
out.snake_Angle = n;
out.snake_Angle_sort = st1;

 [~, g1] = max(n2');
 [~,st1] = sort(g1);
figure(); imagesc((n2(st1,:)));
title('distance snake');
xticks([1 10  20  30 ]);
xticklabels({'0','1.5','3','4.5'})
xlabel('Distance from goal');


out.snake_Distance = n2;
out.snake_Distance_sort = st1;

% significance test.. 
disp( 'Starting Shuffle...');

for i = 1:numROIs;
ROI = i;

% randomization 
for ii = 1: 100;
    rand_ind2try = randi(size(D1,2),size(ind2try,1),1); % circ shift...
    rand_spikes_angle(ii,:) = D1(rand_ind2try);
    rand_spikes_distance(ii,:) = distance_reward1(rand_ind2try);
    
% Angle:
    C1 = histcounts(rand_spikes_angle(ii,:),'Normalization','probability','Binedges',[-180:30:180]);
  
% Distance:
    D1b = histcounts(rand_spikes_distance(ii,:),'Normalization','probability','Binedges',[distance_thresh:100:10000]);

    % smooth data a bit
    norm_angle = C1./B1;
    norm_angle(isnan(norm_angle)) = 0;
    
    norm_distance = D1b./dB1;
    norm_distance(isnan(norm_distance)) = 0;

    rA1 = smooth(norm_angle,angle_smooth); % angle
    rA1_dist = smooth(norm_distance,distance_smooth); % distance
    
    rnd_max_pk_angle(i,ii) = max(rA1)-min(rA1); clear rA1;
    rnd_max_pk_distance(i,ii) = max(rA1_dist)-min(rA1_dist); clear rA1;
    clear rand_ind2try rand_spikes C1
end

% n(i,:) = mat2gray(smooth(A1./B1,10)); clear A1 B1
end

for i = 1:numROIs
p_a(i) = size(find(rnd_max_pk_angle(i,1:100)>rl_pk_angle(i)),2)./size(rnd_max_pk_angle(i,1:100),2);
p_d(i) = size(find(rnd_max_pk_distance(i,1:100)>rl_pk_angle(i)),2)./size(rnd_max_pk_distance(i,1:100),2);
end

out.p_a = p_a;
out.p_d = p_d;




% Plot the current flight paths:
figure(); plot(F2save(:,1), F2save(:,2)) 
title( 'Current flight paths');

% histo the data
figure(); 
subplot(1,3,1);
hold on;
plot(F2save(:,1), F2save(:,2)) 
for i = 1:2;
plot(reward_loc{i}.M(1),reward_loc{i}.M(2),'o','MarkerSize',50,'color',col{i})
title( 'Current flight paths');
end

subplot(1,3,2);
histogram(D1,'Normalization','probability','BinWidth',5);
title(' angles to reward location 1')
xlabel('angle from reward');
ylabel('rel. proportion');
xlim([-180 180]);

subplot(1,3,3);
histogram(D2,'Normalization','probability','BinWidth',5);
xlabel('angle from reward');
ylabel('rel. proportion');
xlim([-180 180]);
title(' angles to reward location 2')


