
function [flightPaths] = ImBat_flightsAngelo_MCS(AllFlights,AllFlightsTime,MicrophoneVect,MicrophoneTime,EcholocationVect,rewind_value,outliers,cluster_to_plot,audioConCat,ROI_Data,varargin)

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

%% Clustering flights
%Cut out flights, downsample to ds_clus positions per flight
all_flights = NaN(3,max(flight_ends-flight_starts_new),num_flights);    %3D matrix with all flights
all_flights_ds = NaN(3,ds_clus,num_flights);                        %3D matrix with all flights(downsampled)
%all_mic_data = NaN(1,max(flight_ends-flight_starts_new),num_flights);

for nf = 1 : size(all_flights,3)
    trajectory = x_spl(:,flight_starts_new(nf):flight_ends(nf));
    velocity = v_abs(:,flight_starts_new(nf):flight_ends(nf));
%    mic_data = MicrophoneVect(flight_starts_new(nf):flight_ends(nf));
    all_flights(:,1:(flight_ends(nf)-flight_starts_new(nf))+1,nf) = trajectory;
    all_flights_vel(1,1:(flight_ends(nf)-flight_starts_new(nf)+1),nf) = velocity;
    all_flights_ds(:,:,nf) = interp1(linspace(1,3,size(trajectory,2)),trajectory',linspace(1,3,ds_clus),'spline')';
 %   all_mic_data(:,1:(flight_ends(nf)-flight_starts_new(nf))+1,nf) = mic_data;
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
    %mic_data_og = MicrophoneVect(flight_starts_new(nf):flight_ends(nf));
    
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
%flight.mic = all_mic_data;
    
% Sort structure according to cluster id
clear flight_sorted;
[flight_sorted.id,I] = sort(flight.id);
flight_sorted.strt_frame = flight.strt_frame(I);
flight_sorted.stop_frame = flight.stop_frame(I);
flight_sorted.pos = flight.pos(:,:,I);
flight_sorted.vel = flight.vel(:,:,I);
flight_sorted.Fs = flight.Fs;
flight_sorted.N = size(flight_sorted.id,1);
%flight_sorted.mic = flight.mic(:,:,I);

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
%flightPaths.mic = flight_sorted.mic;

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
    
    plot3(x_spl(1,:),x_spl(2,:),x_spl(3,:),':','Color',[0.7 0.7 0.7],'MarkerSize',0.001);
    xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
    xlabel('x');    ylabel('y');    zlabel('z');    view(0,0);
    hold off;
    
    subplot(3,n_surv_clusters,n_surv_clusters*2+jj);
    histogram(flightPaths.dur(id));
    xlim([0 15]);   xlabel('Duration(s)');  ylabel('Counts');
end

flightPaths.flightPathsClusterEach = plotFlightPathsClusterEach;

%% Plot the mic and echolocation vectors and put lines where the
% flights are. 
flightpath_starts_binary = NaN(size(flightPaths.trajectoriesSpline,2),1);
flightpath_starts_binary(flightPaths.flight_starts_idx) = .45;
flightpath_ends_binary = NaN(size(flightPaths.trajectoriesSpline,2),1);
flightpath_ends_binary(flightPaths.flight_ends_idx) = .4;
figure(); hold on;
plot(et,MicrophoneVect);
plot(t,flightpath_starts_binary,'og');
plot(t,flightpath_ends_binary,'om');
plot(et,EcholocationVect(:,1),'*r'); plot(et,EcholocationVect(:,2),'*r'); plot(et,EcholocationVect(:,3),'*r'); plot(et,EcholocationVect(:,4),'*r'); plot(et,EcholocationVect(:,5),'*r');

% Pad the EcholocationVect and et vectors to be the size of the
% flightpath vectors. Then segment the echolocation data.
% Collapse EcholocationVect into one
if size(EcholocationVect,2)>1
    aa = EcholocationVect(:,1); 
    bb = EcholocationVect(:,2); 
    cc = EcholocationVect(:,3);
    dd = EcholocationVect(:,4);
    gg = EcholocationVect(:,5);
    aa(isnan(aa)) = bb(isnan(aa)); aa(isnan(aa)) = cc(isnan(aa)); aa(isnan(aa)) = dd(isnan(aa)); aa(isnan(aa)) = gg(isnan(aa));
end

if size(EcholocationVect,2)>1
    Full_EcholocationVect = squeeze(aa);
else
    Full_EcholocationVect = EcholocationVect;
end
   
goal_size = length(t);
current_size = length(et);
% In this case we will just add NaN values to the beginning of the vector
% because we have no echolocation data about that.
nans_to_add = NaN(goal_size-current_size,1);
% This adds this many seconds of nans to the beginning
size(nans_to_add,1)/120;
et_padded = [nans_to_add;et];
EcholocationVect_padded = [nans_to_add;Full_EcholocationVect];
% Now et_padded is the size of goal_size and we can use the flight_starts and flight_ends indexes 

% Replot to make sure NaN values didn't mess it up
flightpath_starts_binary = NaN(size(flightPaths.trajectoriesSpline,2),1);
flightpath_starts_binary(flightPaths.flight_starts_idx) = .45;
flightpath_ends_binary = NaN(size(flightPaths.trajectoriesSpline,2),1);
flightpath_ends_binary(flightPaths.flight_ends_idx) = .4;
figure(); hold on;
plot(et,MicrophoneVect);
plot(t,flightpath_starts_binary,'og');
plot(t,flightpath_ends_binary,'om');
plot(t,EcholocationVect_padded,'*r');
plot(et,Full_EcholocationVect,'*b');

% Do this for audioConCat
nans_to_add = NaN(length(nans_to_add)*1600,1);
audioConCat_padded = [nans_to_add;audioConCat];

% Sort the index flights
flightPaths.flight_starts_idx_sorted = sort(flightPaths.flight_starts_idx);
flightPaths.flight_ends_idx_sorted = sort(flightPaths.flight_ends_idx);

%% Try to plot one flight 
flightnumber_to_plot = 2;

figure(); hold on;
plot3(flightPaths.pos(1,:,flightnumber_to_plot),flightPaths.pos(2,:,flightnumber_to_plot),flightPaths.pos(3,:,flightnumber_to_plot));
plot3(flightPaths.trajectoriesSpline(1,flightPaths.flight_starts_idx_sorted(flightnumber_to_plot):flightPaths.flight_ends_idx_sorted(flightnumber_to_plot)),flightPaths.trajectoriesSpline(2,flightPaths.flight_starts_idx_sorted(flightnumber_to_plot):flightPaths.flight_ends_idx_sorted(flightnumber_to_plot)),flightPaths.trajectoriesSpline(3,flightPaths.flight_starts_idx_sorted(flightnumber_to_plot):flightPaths.flight_ends_idx_sorted(flightnumber_to_plot)),':');
fpe = EcholocationVect_padded(flightPaths.flight_starts_idx_sorted(flightnumber_to_plot):flightPaths.flight_ends_idx_sorted(flightnumber_to_plot));
fpe(fpe==0)=NaN;
fpe_idxs = [flightPaths.flight_starts_idx_sorted(flightnumber_to_plot):flightPaths.flight_ends_idx_sorted(flightnumber_to_plot)];
for i=1:size(fpe,1)
    if fpe(i)==0.5000
        plot3(flightPaths.trajectoriesSpline(1,fpe_idxs(i)),flightPaths.trajectoriesSpline(2,fpe_idxs(i)),flightPaths.trajectoriesSpline(3,fpe_idxs(i)),'*g');
    end
end     
title(strcat("Plotting flight #",num2str(flightnumber_to_plot)));

%% Plot a cluster 
% Plot the specified cluster type with the echolocations over laid
plot_counter=0;
cluster_to_plot = 1;%cluster_to_plot;
figure(); hold on; title(strcat("Echolocations on Cluster ",num2str(cluster_to_plot)," flights"));
flight_subset = find(flightPaths.id==cluster_to_plot);
for i=1:size(flightPaths.id,1)
    if flightPaths.id(i)==cluster_to_plot
        timeline_idx = find(flightPaths.flight_starts_idx_sorted==flightPaths.flight_starts_idx(i));
        plot_counter = plot_counter+1;
    	subplot(round(sqrt(size(flight_subset,1)))+1,round(sqrt(size(flight_subset,1))),plot_counter);
        hold on; 
        title(num2str(timeline_idx));
        plot3(flightPaths.pos(1,:,i),flightPaths.pos(2,:,i),flightPaths.pos(3,:,i));
        plot3(flightPaths.trajectoriesSpline(1,flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i)),flightPaths.trajectoriesSpline(2,flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i)),flightPaths.trajectoriesSpline(3,flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i)),':');
        fpe = EcholocationVect_padded(flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i));
        fpe(fpe==0)=NaN;
        fpe_idxs = [flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i)];
        for j=1:size(fpe,1)
            if fpe(j)==0.5000
                plot3(flightPaths.trajectoriesSpline(1,fpe_idxs(j)),flightPaths.trajectoriesSpline(2,fpe_idxs(j)),flightPaths.trajectoriesSpline(3,fpe_idxs(j)),'*r');
            end
        end  
    end
end

flightPaths.EcholocationVect_padded = EcholocationVect_padded;

save('flightPaths.mat','flightPaths');

%%
% Do a test; play the sound of that flight, plot that flight and
% echolocations
% for i=3:30
%     flightnumber_to_plot = i;
%     figure(); subplot(2,1,1); hold on;
%     plot3(flightPaths.pos(1,:,flightnumber_to_plot),flightPaths.pos(2,:,flightnumber_to_plot),flightPaths.pos(3,:,flightnumber_to_plot));
%     plot3(flightPaths.trajectoriesSpline(1,flightPaths.flight_starts_idx(flightnumber_to_plot):flightPaths.flight_ends_idx(flightnumber_to_plot)),flightPaths.trajectoriesSpline(2,flightPaths.flight_starts_idx(flightnumber_to_plot):flightPaths.flight_ends_idx(flightnumber_to_plot)),flightPaths.trajectoriesSpline(3,flightPaths.flight_starts_idx(flightnumber_to_plot):flightPaths.flight_ends_idx(flightnumber_to_plot)),':');
%     fpe = EcholocationVect_padded(flightPaths.flight_starts_idx(flightnumber_to_plot):flightPaths.flight_ends_idx(flightnumber_to_plot));
%     fpe(fpe==0)=NaN;
%     fpe_idxs = [flightPaths.flight_starts_idx(flightnumber_to_plot):flightPaths.flight_ends_idx(flightnumber_to_plot)];
%     for i=1:size(fpe,1)
%         if fpe(i)==0.5000
%             plot3(flightPaths.trajectoriesSpline(1,fpe_idxs(i)),flightPaths.trajectoriesSpline(2,fpe_idxs(i)),flightPaths.trajectoriesSpline(3,fpe_idxs(i)),'*g');
%         end
%     end     
%     subplot(2,1,2);
%     plot(flightPaths.vel(1,:,flightnumber_to_plot));
%     title(strcat("Plotting flight #",num2str(flightnumber_to_plot)));
% 
%     sound(audioConCat_padded(flightPaths.flight_starts_idx(flightnumber_to_plot)*1600:flightPaths.flight_ends_idx(flightnumber_to_plot)*1600),fs);
%     pause;
% end
% 
% %% With sorted flight indexes
% ff=100;
% for i=1:30
%     flightnumber_to_plot = i;
%     figure(); subplot(3,1,1); hold on;
%     %plot3(flightPaths.pos(1,:,flightnumber_to_plot),flightPaths.pos(2,:,flightnumber_to_plot),flightPaths.pos(3,:,flightnumber_to_plot));
%     plot3(flightPaths.trajectoriesSpline(1,flightPaths.flight_starts_idx_sorted(flightnumber_to_plot)-ff:flightPaths.flight_ends_idx_sorted(flightnumber_to_plot)+ff),flightPaths.trajectoriesSpline(2,flightPaths.flight_starts_idx_sorted(flightnumber_to_plot)-ff:flightPaths.flight_ends_idx_sorted(flightnumber_to_plot)+ff),flightPaths.trajectoriesSpline(3,flightPaths.flight_starts_idx_sorted(flightnumber_to_plot)-ff:flightPaths.flight_ends_idx_sorted(flightnumber_to_plot)+ff),':');
%     fpe = EcholocationVect_padded(flightPaths.flight_starts_idx_sorted(flightnumber_to_plot)-ff:flightPaths.flight_ends_idx_sorted(flightnumber_to_plot)+ff);
%     fpe(fpe==0)=NaN;
%     fpe_idxs = [flightPaths.flight_starts_idx_sorted(flightnumber_to_plot)-ff:flightPaths.flight_ends_idx_sorted(flightnumber_to_plot)+ff];
%     for j=1:size(fpe,1)
%         if fpe(j)==0.5000
%             plot3(flightPaths.trajectoriesSpline(1,fpe_idxs(j)),flightPaths.trajectoriesSpline(2,fpe_idxs(j)),flightPaths.trajectoriesSpline(3,fpe_idxs(j)),'*g');
%         end
%     end     
%     title(strcat("Plotting flight #",num2str(flightnumber_to_plot)));
% 
%     subplot(3,1,2); hold on; 
%     plot(audioConCat_padded(flightPaths.flight_starts_idx_sorted(flightnumber_to_plot)*1600:flightPaths.flight_ends_idx_sorted(flightnumber_to_plot)*1600));
%     subplot(3,1,3); hold on;
%     plot(flightPaths.batSpeed(flightPaths.flight_starts_idx_sorted(flightnumber_to_plot):flightPaths.flight_ends_idx_sorted(flightnumber_to_plot)));
%     sound(audioConCat_padded(flightPaths.flight_starts_idx_sorted(flightnumber_to_plot)*1600:flightPaths.flight_ends_idx_sorted(flightnumber_to_plot)*1600),fs);
%     pause;
% end
end
