
function [flightPaths] = ImBat_MCS_flightsAngelo_offset(AllFlights,AllFlightsTime,MicrophoneVect,MicrophoneTime,EcholocationVect,rewind_value,outliers,cluster_to_plot,audioConCat,ROI_Data,varargin)

% 

% Default paramaters:

% Clustering params
ds_clus = 6;                                                                %number of 3D-points/flight for clustering 
%madeleine 25 splines, PCA+, 1m linkage
%angelo 6 splines, PCA-, 0.7m linakge, min 5 
pca_features = false;                                                       %if using PCA
k_means = false;                                                            %if using k-means
dist = 1.5;                                                                 %linkage distance
reassign = true;                                                            %re-order clusters according to density
N_min = 3;  
DayIndex = ones(size(AllFlightsTime,1),1);                                    % if only clustering on one day
VideoFrameRate = 120;
%EcholocationVect = EcholocationVect';

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'fs'
            fs=varargin{i+1};
        case 'dist'
            dist=varargin{i+1};
        case 'n_min'
            N_min=varargin{i+1};
        case 'n_splines'
            ds_clus = varargin{i+1};
        case 'day_index'
            day_index = varargin{i+1};
        case 'joint'
            joint = varargin{i+1};
    end
end

% General Params
CNMF_Fs = fs;                                                              % Acquisition frequency (Hz) Ca-imaging (after CNMF-e)
t = AllFlightsTime;
et = MicrophoneTime;

% % Flight Room references (provvisory)
xR = +2.85; xL = -2.85; yF = 2.50;  yB = -2.50;  zT = 2.20;                 %Flight volume coordinates
F3 = [2.56; 1.43; 0.72];    F4 = [2.56; -1.24; 0.72];                       %Feeders coordinates
% 
% T_tk = size(trackData.Markers,1);  t = [0:1/trackData.VideoFrameRate:(T_tk-1)/trackData.VideoFrameRate];  %Generate time vector
% rew_signal = trackData.AnalogSignals(:,1);                                            %Reward signal, usually on Analog Ch1


%% Trajectory extraction
x_mean = [ AllFlights(:,1)  AllFlights(:,2)  AllFlights(:,3)]'./1000;  

% Filter and interpolate
x_filt = x_mean; %medfilt1(x_mean,VideoFrameRate/2,[],2,'omitnan','truncate'); %filter after interpolating
x_intr = x_mean;% fillmissing(x_filt,'next',2,'EndValues','nearest');
x_spl = x_mean; %x_intr; %csaps(t, x_intr, 0.9, t);

%Frame rate
new_t = t';
tracking_Fs = VideoFrameRate;

%% Calculate velocity and flight segmentation
v = diff(x_spl,1,2)./[diff(new_t);diff(new_t);diff(new_t)]; v=[zeros(3,1) v];
v_abs = vecnorm(v,2,1);

nonflying = find(v_abs < 1);        
toofast = find(v_abs > 30);
x_flying = x_spl;                   x_flying(:,[nonflying toofast]) = nan;
batspeed = v_abs;                   batspeed([nonflying toofast]) = nan;
bflying=~isnan(batspeed)';           %vector of 1s when the bat is flying

% For each sample, sum up the next 1s of data(flights are longer than 1s),Code from Nick D.
allsums = [];
for bf = 1 : size(bflying,1)-tracking_Fs
    allsums(bf) = sum(bflying(bf:bf+tracking_Fs));
end

% Detect flight starts and stops
[R,rLT,rUT,rLL,rUL] = risetime(allsums);    
[F,fLT,fUT,fLL,fUL] = falltime(allsums);           
if length(R) ~= length(F)
fLT(length(R)) = length(allsums);
fUT(length(R)) = length(allsums);
F(length(R)) = F(length(F));
end
flight_starts = round(rLT+tracking_Fs/2);
flight_starts_new = round(rLT+tracking_Fs/2);%-rewind_value;
flight_ends = round(fLT+tracking_Fs/2); %... +Fs is a sistematic correction, useful
num_flights = size(R,2);
ref = ones(size(flight_starts));
avg_flight_time = mean((flight_ends-flight_starts)./tracking_Fs);

% Plot the velocity at the normal flight starts and at the rewound flight starts


%% Plot results

% Plot 2D flight trajectories
plotFlightPathsAll = figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
plot(x_mean(1,:),x_mean(2,:),'.');
hold on;        rectangle('Position',[xL yB xR-xL yF-yB]);
scatter([F3(1) F4(1)],[F3(2) F4(2)],'filled');  hold off;
xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
title(['Raw flights']);
xlabel('m'); ylabel('m');
hold off

subplot(1,2,2);
plot(x_spl(1,:),x_spl(2,:)); hold on; plot(x_mean(1,:),x_mean(2,:),'.','MarkerSize',1);
rectangle('Position',[xL yB xR-xL yF-yB]);  hold off;
xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
title(['Spline flights: ']);
xlabel('m'); ylabel('m');
hold off

% Plot session timeline
plotFlightTimeline = figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
ax1 = subplot(3,1,1);   plot(t,x_mean(1,:),'.');  hold on;
plot(new_t,x_spl(1,:),'--');   refline(0,F3(1));    hold off;
legend('cluster/mean','spl');   ylabel('x (m)');
ax2 = subplot(3,1,2);   plot(new_t,v_abs,'.');
hold on;    stem(new_t(flight_starts),ref);    stem(new_t(flight_ends),ref);  hold off;
ylabel('v (m/s)');
 ax3 = subplot(3,1,3);  % plot(t,rew_signal);
ylabel('Rewards');
linkaxes([ax1,ax2,ax3],'x');    xlabel('Samples');

% MCS Plot session timeline
[new_t_s,new_t_s_idx] = sort(new_t);
figure(); hold on; plot(new_t(new_t_s_idx),v_abs,'.'); 
fs_new_t = new_t(new_t_s_idx);
stem(fs_new_t(flight_starts),ref); 
disp(fs_new_t(flight_starts(1:5)));

% MCS Plot Zoomed in outlier velocity traces
flight_starts_new = flight_starts;
for out_i=1:size(outliers,2)
    i=outliers(out_i);
   
    %figure(); hold on; title(strcat("Outlier ",num2str(out_i))); 
    st_temp = fs_new_t(flight_starts(i));
    %plot((fs_new_t(flight_starts(i)-rewind_value*2:flight_starts(i+1))),v_abs(flight_starts(i)-rewind_value*2:flight_starts(i+1)),'.'); 
    %stem(fs_new_t(flight_starts(i)),1);
    %stem(fs_new_t(flight_starts(i)-rewind_value),1);
    new_outlier_flight_start = fs_new_t(flight_starts(i)-rewind_value);
    for j=1:max(flight_starts)
        if fs_new_t(j) == new_outlier_flight_start
            new_temp_flight_start = j;
            flight_starts_new(outliers(out_i)) = new_temp_flight_start;
        end
    end
    if i == size(flight_starts,2)
        continue
    else
        figure(); hold on; title(strcat("New Flight Start ",num2str(out_i))); 
        st_temp = fs_new_t(flight_starts(i));
        plot((fs_new_t(flight_starts(i)-rewind_value*2:flight_starts(i+1))),v_abs(flight_starts(i)-rewind_value*2:flight_starts(i+1)),'.'); 
        stem(fs_new_t(flight_starts_new(i)),1);
    end
end

%flight_starts = flight_starts_new;
% Plot flights in color time order
plotFlightPathsStartStop = figure();
if size(R,2) > 0
    CM = jet(size(R,2));
    for nf = 1 : size(R,2)
        hold on
        plot3(x_spl(1,flight_starts_new(nf):flight_ends(nf)),x_spl(2,flight_starts_new(nf):flight_ends(nf)),x_spl(3,flight_starts_new(nf):flight_ends(nf)),'LineWidth',1,'Color',CM(nf,:))
        hold on
        
        fstartxyz(nf,1) = x_spl(1,flight_starts_new(nf)); %round(nanmean(x_spl(1,flight_starts(nf):flight_starts(nf)+trackData.VideoFrameRate/2)));
        fstartxyz(nf,2) = x_spl(2,flight_starts_new(nf)); %round(nanmean(x_spl(2,flight_starts(nf):flight_starts(nf)+trackData.VideoFrameRate/2)));
        fstartxyz(nf,3) = x_spl(3,flight_starts_new(nf)); %round(nanmean(x_spl(3,flight_starts(nf):flight_starts(nf)+trackData.VideoFrameRate/2)));
        
        fendxyz(nf,1) = x_spl(1,flight_ends(nf)); %round(nanmean(x_spl(1,flight_ends(nf):flight_ends(nf)+trackData.VideoFrameRate/2)));
        fendxyz(nf,2) = x_spl(2,flight_ends(nf)); %round(nanmean(x_spl(2,flight_ends(nf):flight_ends(nf)+trackData.VideoFrameRate/2)));
        fendxyz(nf,3) = x_spl(3,flight_ends(nf)); %round(nanmean(x_spl(3,flight_ends(nf):flight_ends(nf)+trackData.VideoFrameRate/2)));
        
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
title(['All flights start(r)/stop(b): ']);
% modify labels for tick marks
view(0,90)
xlim([-3 3])
ylim([-3 3])
xlabel('m'); ylabel('m');

hold off

% Plot an example flight with the microphone trace underneath it
% Can't define MicrophoneVect or EcholocationVect with the indexes of
% flight_starts (different vector lenghts)
figure();
subplot(2,1,1); hold on;
plot3(x_spl(1,flight_starts_new(1):flight_ends(1)),x_spl(2,flight_starts_new(1):flight_ends(1)),x_spl(3,flight_starts_new(1):flight_ends(1)),'LineWidth',1,'Color',CM(1,:))
subplot(2,1,2);  hold on;
plot(MicrophoneVect(flight_starts_new(1):flight_ends(1)));
plot(EcholocationVect(flight_starts_new(1):flight_ends(1)),'*r');

% Plot each flight the with echolocation pips on the trace
% One example:
figure();
hold on;
plot3(x_spl(1,flight_starts_new(1):flight_ends(1)),x_spl(2,flight_starts_new(1):flight_ends(1)),x_spl(3,flight_starts_new(1):flight_ends(1)),'LineWidth',1,'Color',CM(1,:))
echo_idxs = find(EcholocationVect==0.01);
echo_1 = echo_idxs(echo_idxs>flight_starts_new(1));
echo_2 = echo_1(echo_1<flight_ends(1));
scatter3(x_spl(1,echo_2),x_spl(2,echo_2),x_spl(3,echo_2),'*r');
scatter3(x_spl(1,flight_starts_new(1)),x_spl(2,flight_starts_new(1)),x_spl(3,flight_starts_new(1)),'ob');


%% Clustering flights
%Cut out flights, downsample to ds_clus positions per flight
all_flights = NaN(3,max(flight_ends-flight_starts_new),num_flights);    %3D matrix with all flights
all_flights_ds = NaN(3,ds_clus,num_flights);                        %3D matrix with all flights(downsampled)
all_echolocations =  NaN(1,max(flight_ends-flight_starts_new),num_flights);
all_mic_data = NaN(1,max(flight_ends-flight_starts_new),num_flights);

for nf = 1 : size(all_flights,3)
    trajectory = x_spl(:,flight_starts_new(nf):flight_ends(nf));
    velocity = v_abs(:,flight_starts_new(nf):flight_ends(nf));
    echolocations = EcholocationVect(flight_starts_new(nf):flight_ends(nf));
    mic_data = MicrophoneVect(flight_starts_new(nf):flight_ends(nf));
    try
        joint_echolocations = JT(:,flight_starts_new(nf):flight_ends(nf));
    catch
        disp("No Joint Echolocation Vector loaded");
    end
    all_flights(:,1:(flight_ends(nf)-flight_starts_new(nf))+1,nf) = trajectory;
    all_flights_vel(1,1:(flight_ends(nf)-flight_starts_new(nf)+1),nf) = velocity;
    all_flights_ds(:,:,nf) = interp1(linspace(1,3,size(trajectory,2)),trajectory',linspace(1,3,ds_clus),'spline')';
    all_echolocations(:,1:(flight_ends(nf)-flight_starts_new(nf))+1,nf) = echolocations;
    all_mic_data(:,1:(flight_ends(nf)-flight_starts_new(nf))+1,nf) = mic_data;
    try
        all_joint_echolocations(:,1:(flight_ends(nf)-flight_starts_new(nf))+1,nf) = joint_echolocations;
    catch
        disp("No Joint Echolocation Vector loaded");
    end
    % MCS added
    trajectory_og = x_spl(:,flight_starts(nf):flight_ends(nf));
    velocity_og = v_abs(:,flight_starts(nf):flight_ends(nf));
    all_flights_og(:,1:(flight_ends(nf)-flight_starts(nf))+1,nf) = trajectory_og;
    all_flights_vel_og(1,1:(flight_ends(nf)-flight_starts(nf)+1),nf) = velocity_og;
    all_flights_ds_og(:,:,nf) = interp1(linspace(1,3,size(trajectory_og,2)),trajectory_og',linspace(1,3,ds_clus),'spline')';
    mic_data_og = MicrophoneVect(flight_starts_new(nf):flight_ends(nf));
    
    %Uncomment if you want to see how the downsampled flights look like
%     figure(); hold on;
%     plot3(all_flights(1,:,nf),all_flights(2,:,nf),all_flights(3,:,nf),'Color','b');
%     plot3(all_flights_ds(1,:,nf),all_flights_ds(2,:,nf),all_flights_ds(3,:,nf),'Color','r');
%     
%     plot3(all_flights_og(1,:,nf),all_flights_og(2,:,nf),all_flights_og(3,:,nf),'Color','g');
%     plot3(all_flights_ds_og(1,:,nf),all_flights_ds_og(2,:,nf),all_flights_ds_og(3,:,nf),'Color','c');
    
end

% Define X matrix of features for clustering (downsampled coordinates, stacked together)
X = [all_flights_ds(1,:,:), all_flights_ds(2,:,:), all_flights_ds(3,:,:)];
X = reshape(X,3*size(all_flights_ds,2),size(R,2));
X = X';     %so then X = #flights x #features

% If dimensionality reduction is needed
if pca_features
    [coeff,score,latent] = pca(X);     X = score(:,1:5);
end

% k-means or hierarchical clustering (with euclidean distance and shortest linkage)
if k_means
    n_clusters = 15;    idx = kmeans(X,n_clusters);
else
    plotClusterDistance = figure();
    Y = pdist(X,'euclidean');   Z = linkage(Y);
    hLines = dendrogram(Z,0);  hold on;    refline(0,dist);    hold off;
    idx = cluster(Z,'Cutoff',dist,'Criterion','distance');
    title([num2str(length(unique(idx))) ' clusters: ']);
    ylim([0 10]);
end

% Create structure with flight start stop frames, id of the trajectory
clear flight 
flight.strt_frame = ceil(flight_starts_new)';
flight.stop_frame = ceil(flight_ends)';
flight.pos = all_flights;
flight.vel = all_flights_vel;
flight.id = idx;
flight.Fs = tracking_Fs;
flight.echo = all_echolocations;
flight.mic = all_mic_data;
try
    flight.JT = all_joint_echolocations;
catch
    disp("No JT");
end

no_echos=0;
for i=1:size(flight.echo,3)
    if nansum(all_echolocations(1,:,i))*100 == 0
        no_echos = no_echos+1;
    end
end
disp(strcat("Considering one microphone, ",num2str(100-(no_echos/size(flight.echo,3)*100)),"% of flights have echolocations"));

try
    no_echos=0;
    for i=1:size(flight.JT,3)
        if nansum(all_joint_echolocations(1,:,i))*100 == 0
            no_echos = no_echos+1;
        end
    end
    disp(strcat("Considering ALL mics, ",num2str(100-(no_echos/size(flight.JT,3)*100)),"% of flights have echolocations"));
catch
    disp("Only one mic of data");
end
    

% Add in the echolocation pip vectors
% echo_frames = NaN(1,50);
% echo_idxs = find(EcholocationVect==0.01);
% for i=1:size(flight_starts_new,2)
%     echo_1 = echo_idxs(echo_idxs>flight_starts_new(i));
%     echo_2 = echo_1(echo_1<flight_ends(i));
%     echo_2_padded = [echo_2,NaN(1,50-size(echo_2,2))];
%     echo_frames(i,:) = echo_2_padded;
% end
% flight.echo_frames = echo_frames;

% Sort structure according to cluster id
clear flight_sorted;
[flight_sorted.id,I] = sort(flight.id);
flight_sorted.strt_frame = flight.strt_frame(I);
flight_sorted.stop_frame = flight.stop_frame(I);
flight_sorted.pos = flight.pos(:,:,I);
flight_sorted.vel = flight.vel(:,:,I);
flight_sorted.Fs = flight.Fs;
flight_sorted.N = size(flight_sorted.id,1);
flight_sorted.echos = flight.echo(:,:,I);
flight_sorted.mic = flight.mic(:,:,I);
try
    flight_sorted.JT = flight.JT(:,:,I);
catch
    disp("No JT");
end

% Assign isolated clusters to cluster #flights+1
[Ns,b] = histc(flight_sorted.id,unique(flight_sorted.id));
flight_sorted.id(Ns(b)<N_min) = size(all_flights,3)+1;             %flight_sorted.id(Ns(b)==1) = size(all_flights,3)+1;
id_surv_clusters = unique(flight_sorted.id);
n_surv_clusters = size(id_surv_clusters,1);

% Create final structure flight.clus after re-assignment
clear flightPaths;
flightPaths.id = flight_sorted.id;
flightPaths.flight_starts_idx = flight_sorted.strt_frame';
flightPaths.flight_ends_idx = flight_sorted.stop_frame';
flightPaths.pos = flight_sorted.pos;
flightPaths.vel = flight_sorted.vel;
flightPaths.Fs = flight_sorted.Fs;
flightPaths.N = flight_sorted.N;
flightPaths.day = day_index(flightPaths.flight_starts_idx);% this is the day index...
flightPaths.echos = flight_sorted.echos;
flightPaths.mic = flight_sorted.mic;
try
    flightPaths.JT = flight_sorted.JT;
catch
    disp("No JT");
end

for jj=1:n_surv_clusters;
    flightPaths.id(flight_sorted.id == id_surv_clusters(jj)) = jj;
end
id_surv_clusters = unique(flightPaths.id);

%Re-assign id for convenience, if necessary
if reassign
    new_ord = [];
    [~,new_ord] = sort(histc(flightPaths.id,id_surv_clusters(1:end-1)),'descend');
    new_ord = [new_ord; id_surv_clusters(end)];
    new_ord = circshift(new_ord,1);
    reassign_matrix =(flightPaths.id == new_ord');
    for jj=1:n_surv_clusters;
        flightPaths.id(reassign_matrix(:,jj)) = jj;
    end  
end 

%Calculate trajectory length, duration in s and interflight (take-off to take-off)
for ii = 1:flight_sorted.N
    flightPaths.length(ii)= arclength(flightPaths.pos(1,~isnan(flightPaths.pos(1,:,ii)),ii),flightPaths.pos(2,~isnan(flightPaths.pos(2,:,ii)),ii),flightPaths.pos(3,~isnan(flightPaths.pos(3,:,ii)),ii));
    flightPaths.dur(ii) = (flightPaths.flight_ends_idx(ii)-flightPaths.flight_starts_idx(ii))./flightPaths.Fs;
end
flightPaths.ifd = diff(flightPaths.flight_starts_idx)';


%group the cluster ids
for i = 1:max(flightPaths.id)
    flightPaths.clusterIndex{i} = find(flightPaths.id == i);
end

%add Tobias specific variables for future plots
flightPaths.trajectoriesContinous = x_intr;
flightPaths.trajectoriesSpline = x_spl;
flightPaths.tracjectoriesRaw = x_mean;
flightPaths.batSpeed = v_abs';
flightPaths.flight_starts_xyz = flightPaths.pos(:,1,:); %starting position of each flight
flightPaths.flight_ends_xyz = zeros(length(flightPaths.pos(1,1,:)),3); %make matrix for landing position for each flight 
for i = 1:length(flightPaths.pos(1,1,:))
    try
        flightPaths.flight_ends_xyz(i,:) = flightPaths.pos(:,find(isnan(flightPaths.pos(1,:,i)),1)-1,i); %find last xyz position
    catch
        flightPaths.flight_ends_xyz(i,:) = flightPaths.pos(:,end,i);
    end
end
flightPaths.flightTimeline = plotFlightTimeline;
flightPaths.flightPathsAll = plotFlightPathsAll;
flightPaths.flightPathsStartStop = plotFlightPathsStartStop;
flightPaths.clusterDistance = plotClusterDistance;
flightPaths.ds_clus = ds_clus;                                                                %number of 3D-points/flight for clustering 
flightPaths.pca_features = pca_features;                                                       %if using PCA
flightPaths.linkDist = dist;                                                                 %linkage distance
flightPaths.N_min = N_min;
%% Visualize

% Force min cluster: 
n_surv_clusters = 5;

plotFlightPathsClusterEach = figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
col = hsv(n_surv_clusters);
title(['Flight clusters:']);
for jj=1:n_surv_clusters;
    id = find(flightPaths.id==jj);
    
    subplot(3,n_surv_clusters,jj);
    avg_take_off = [];
    for ii=1:size(id,1);
        hold on;
        title(['Cluster' num2str(jj) '  (' num2str(size(id,1)) ' flights)'])
        plot3(flightPaths.pos(1,:,id(ii)),flightPaths.pos(2,:,id(ii)),flightPaths.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col(jj,:));
        avg_take_off = [avg_take_off flightPaths.pos(:,1,id(ii))];
        hold on;
    end
    take_off = mean(avg_take_off,2);
    %textscatter3(take_off(1),take_off(2),take_off(3),"Take-off");
    
    plot3(x_spl(1,:),x_spl(2,:),x_spl(3,:),':','Color',[0.7 0.7 0.7],'MarkerSize',0.001);
    xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
    xlabel('x');    ylabel('y');    zlabel('z');    view(2);
    hold off;
    
    subplot(3,n_surv_clusters,n_surv_clusters+jj);
    avg_take_off = [];
    for ii=1:size(id,1);
        hold on;
        title(['Cluster' num2str(jj) '  (' num2str(size(id,1)) ' flights)'])
        plot3(flightPaths.pos(1,:,id(ii)),flightPaths.pos(2,:,id(ii)),flightPaths.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col(jj,:));
        avg_take_off = [avg_take_off flightPaths.pos(:,1,id(ii))];
        hold on;
    end
    take_off = mean(avg_take_off,2);
    %textscatter3(take_off(1),take_off(2),take_off(3),"Take-off");
    
    plot3(x_spl(1,:),x_spl(2,:),x_spl(3,:),':','Color',[0.7 0.7 0.7],'MarkerSize',0.001);
    xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
    xlabel('x');    ylabel('y');    zlabel('z');    view(0,0);
    hold off;
    
    subplot(3,n_surv_clusters,n_surv_clusters*2+jj);
    histogram(flightPaths.dur(id));
    xlim([0 15]);   xlabel('Duration(s)');  ylabel('Counts');
end

flightPaths.flightPathsClusterEach = plotFlightPathsClusterEach;

% Plot the specified cluster type with the echolocations over laid
% If there is a JT matrix for this (flight_echoVect_joint) use this!
cluster_to_plot = 1%cluster_to_plot;
figure(); title(strcat("Echolocations on Cluster ",num2str(cluster_to_plot)," flights"));
flight_subset = find(flightPaths.id==cluster_to_plot);
for i=1:size(flight_subset,1)
    EF = flightPaths.echos(:,:,flight_subset(i));
    ef = EF(~isnan(EF));
    ef_2 = ef(ef~=0);
    subplot(round(sqrt(size(flight_subset,1)))+1,round(sqrt(size(flight_subset,1))),i);
    hold on; 
    if isempty(ef_2)
        title('No Echolocations for this flight');
    end
    %for j=1:size(ef_2,2)
        fpe = squeeze(flightPaths.echos(:,:,flight_subset(i)));
        fpe(fpe==0)=NaN;
        fpe_loc = find(~isnan(fpe));
%         fidx_1 = flightPaths.flight_starts_idx-ef(j);
%         fidx_1(fidx_1 < 0) = 1000000;
%         fidx_2 = find(fidx_1 == min(fidx_1));
%         fidx_3 = find(flightPaths.flight_starts_idx == flightPaths.flight_starts_idx(fidx_2));
        plot3(flightPaths.pos(1,:,flight_subset(i)),flightPaths.pos(2,:,flight_subset(i)),flightPaths.pos(3,:,flight_subset(i)));
        scatter3(flightPaths.pos(1,fpe_loc,flight_subset(i)), flightPaths.pos(2,fpe_loc,flight_subset(i)),flightPaths.pos(3,fpe_loc,flight_subset(i)),'*r');
    %end
end


% Try JT thing
% cluster_to_plot = cluster_to_plot;
% figure(); title(strcat("Echolocations from all mics on Cluster ",num2str(cluster_to_plot)," flights"));
% flight_subset = find(flightPaths.id==cluster_to_plot);
% for i=1:size(flight_subset,1)
%     EF = flightPaths.JT(:,:,flight_subset(i));
%     ef = EF(~isnan(EF));
%     ef_2 = ef(ef~=0);
%     subplot(round(sqrt(size(flight_subset,1)))+1,round(sqrt(size(flight_subset,1))),i);
%     hold on; 
%     if isempty(ef_2)
%         title('No Echolocations for this flight');
%     end
%     %for j=1:size(ef_2,2)
%         fpe = squeeze(flightPaths.JT(:,:,flight_subset(i)));
%         fpe_loc = find(fpe==0.01);
% %         fidx_1 = flightPaths.flight_starts_idx-ef(j);
% %         fidx_1(fidx_1 < 0) = 1000000;
% %         fidx_2 = find(fidx_1 == min(fidx_1));
% %         fidx_3 = find(flightPaths.flight_starts_idx == flightPaths.flight_starts_idx(fidx_2));
%         plot3(flightPaths.pos(1,:,flight_subset(i)),flightPaths.pos(2,:,flight_subset(i)),flightPaths.pos(3,:,flight_subset(i)));
%         scatter3(flightPaths.pos(1,fpe_loc,flight_subset(i)), flightPaths.pos(2,fpe_loc,flight_subset(i)),flightPaths.pos(3,fpe_loc,flight_subset(i)),'*r');
%     %end
% end

% Test- plot the mic and echolocation vectors and put lines where the
% flights are

% %Upsample echolocation vect
% EV = find(EcholocationVect==0.5);
% EV_UP = EV*1600;
% EE = NaN(size(audioConCat,1),1);
% EE(EV_UP) = 0.5;
% 
% y = upsample(ROI_Data{1}.Alignment.metrics.ttl_tv_offset,1600)
% ttl_offset_zero_point = find(ROI_Data{1}.Alignment.metrics.ttl_tv_offset == 0);
% 
% % Offset things?
% ACC = [audioConCat(ttl_offset_zero_point+1:end);NaN(ttl_offset_zero_point,1)];
% EEE = [EE(ttl_offset_zero_point+1:end);NaN(ttl_offset_zero_point,1)];
% for i=2:size(flightPaths.pos,3) 
%     FSI(i) = (flightPaths.flight_starts_idx(i)-ttl_offset_zero_point)*1600;
%     FSTI(i) = (flightPaths.flight_ends_idx(i)-ttl_offset_zero_point)*1600;
% end
% FSVect = NaN(size(EEE,1),1);
% FSTIVect = NaN(size(EEE,1),1);
% FSVect(FSI) = 1;
% FSTIVect(FSTI) = 1;
% 
% figure(); hold on; 
% plot(ACC(1000000:2000000));%(1:round(size(EE,1)/4)));
% plot(EEE(1000000:2000000),'*r');%(1:round(size(EE,1)/4)),'*r');
% for i=2:size(flightPaths.pos,3) 
%     FSI = (flightPaths.flight_starts_idx(i)-ttl_offset_zero_point)*1600;
%     FSTI = (flightPaths.flight_ends_idx(i)-ttl_offset_zero_point)*1600;
%     xline(FSI,'g','LineWidth',3);
%     xline(FSTI,'m');   
% end

% Turn flightpaths starts into a binary vector
% Plots this

figure(); hold on;
title("All Flight Starts and Stops Aligned with Echolocation Vector");
for i=2:size(flightPaths.pos,3) 
    plot(flightPaths.flight_starts_idx(i)-ttl_offset_zero_point));
%     FSTI = (flightPaths.flight_ends_idx(i)-ttl_offset_zero_point)*1600;


% Do a test; play the sound of that flight, plot that flight and
% echolocations
for i=2:30
flightnumber_to_plot = i;
fpe = squeeze(flightPaths.echos(:,:,flightnumber_to_plot));
fpe(fpe==0)=NaN;
fpe_loc = find(~isnan(fpe));

ee = find(isnan(flightPaths.pos(1,:,flightnumber_to_plot)));
flight_length = ee(1)-1;
flight_length_upsampled = flight_length*1600;
flight_starting_index = 1;
sound(audioConCat(flightPaths.flight_starts_idx(flightnumber_to_plot)*1600:flightPaths.flight_ends_idx(flightnumber_to_plot)*1600),fs);

figure(); hold on;
plot3(flightPaths.pos(1,:,flightnumber_to_plot),flightPaths.pos(2,:,flightnumber_to_plot),flightPaths.pos(3,:,flightnumber_to_plot));
scatter3(flightPaths.pos(1,fpe_loc,flightnumber_to_plot), flightPaths.pos(2,fpe_loc,flightnumber_to_plot),flightPaths.pos(3,fpe_loc,flightnumber_to_plot),'*r');
pause;
sound(audioConCat(flightPaths.flight_ends_idx(flightnumber_to_plot)*1600:flightPaths.flight_starts_idx(flightnumber_to_plot+1)*1600),fs);
pause;
end

%save([analysis_folder '\' batName '_' dateSesh '_' sessionType '_flightPaths.mat'],'flightPaths')
end
