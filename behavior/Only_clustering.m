%% Script for clustering. Run on a folder with a _track.mat file
clear; close;
% General Params
CNMF_Fs = 6;                                                                %Acquisition frequency (Hz) Ca-imaging (after CNMF-e)
analyze_Ca = false;                                                         %If Ca analysis is included
align_to_imaging = false;                                                   %If aligning to imaging
use_bat_cluster = true;                                                     %If using bat_cluster
dwn_sample = false;                                                         %down-sample to CNMF-e
save_data = false;                                                           %save after extraction
save_img = false;                                                            %save cluster image

% Clustering params
ds_clus = 6;                                                                %number of 3D-points/flight for clustering
pca_features = false;                                                       %if using PCA
k_means = false;                                                            %if using k-means
dist = 0.7;                                                                 %linkage distance
reassign = true;                                                            %re-order clusters according to density
N_min = 4;                                                                  %min number of flights for being a cluster

% Flight Room references (provvisory)
xR = +2.85; xL = -2.85; yF = 2.50;  yB = -2.50;  zT = 2.20;                 %Flight volume coordinates
F3 = [2.56; 1.43; 0.72];    F4 = [2.56; -1.24; 0.72];                       %Feeders coordinates

% Create folder for storing the results
analysis_directory=fullfile(pwd,['Analysis_',datestr(now, 'yymmdd_HHMM')]);
if ~exist(analysis_directory,'dir')
    mkdir(analysis_directory);
end

% Load data after CNMF-e extraction
if analyze_Ca
    extracted_ROIs_file = dir(fullfile(cd, 'CNMFE*'));          load(extracted_ROIs_file.name);
end

% Load data from filename_track.mat and align, if needed
tracking_file = dir(fullfile(cd, '*_track.mat'));           load(tracking_file.name); batdate = tracking_file.name(1:10);

T_tk = size(Markers,1);  t = [0:1/VideoFrameRate:(T_tk-1)/VideoFrameRate];  %Generate time vector
rew_signal = AnalogSignals(:,1);                                            %Reward signal, usually on Analog Ch1

%% Trajectory extraction
if use_bat_cluster
    x_mean = [Markers(:,1,1) Markers(:,1,2) Markers(:,1,3)]'./1000;     x_mean(x_mean==0) = nan;
else
    Markers_nan = Markers;  Markers_nan(Markers_nan==0) = nan;
    Markers_mean = mean(Markers_nan,2,'omitnan');
    x_mean = [Markers_mean(:,1,1) Markers_mean(:,1,2) Markers_mean(:,1,3)]'./1000;
end

% Filter and interpolate
x_filt = medfilt1(x_mean,VideoFrameRate/2,[],2,'omitnan','truncate');
x_intr = fillmissing(x_filt,'next',2,'EndValues','nearest');
x_spl = csaps(t, x_intr, 0.9, t);

%% Temporal alignment
if align_to_imaging
    alignment_file = dir(fullfile(cd, 'Alignment*'));      load(alignment_file.name);
    c3d_times = out.Location_time';                    img_times = out.video_times';
    c3d_sampling = diff(c3d_times);                    img_sampling = diff(img_times);
    c3d_sampl_interval = mean(c3d_sampling);           img_sampl_interval = mean(img_sampling);
    c3d_Fs = round(1/c3d_sampl_interval);              img_Fs = round(1/img_sampl_interval);
    
    t = c3d_times;
    new_t = img_times;
    tracking_Fs = img_Fs;
    x_spl = csaps(t, x_spl, 0.98, new_t);
else
    new_t = t;
    tracking_Fs = VideoFrameRate;
end

%% Calculate velocity and flight segmentation
v = diff(x_spl,1,2)./[diff(new_t);diff(new_t);diff(new_t)]; v=[zeros(3,1) v];
v_abs = vecnorm(v,2,1);

nonflying = find(v_abs < 1);        toofast = find(v_abs > 30);
x_flying = x_spl;                   x_flying(:,[nonflying toofast]) = nan;
batspeed = v_abs;                   batspeed([nonflying; toofast]) = nan;
bflying=~isnan(batspeed)';           %vector of 1s when the bat is flying

% For each sample, sum up the next 1s of data(flights are longer than 1s),Code from Nick D.
allsums = [];
for bf = 1 : size(bflying,1)-tracking_Fs
    allsums(bf) = sum(bflying(bf:bf+tracking_Fs));
end

% Detect flight starts and stops
[R,rLT,rUT,rLL,rUL] = risetime(allsums);    flight_starts = round(rLT+tracking_Fs/2);
[F,fLT,fUT,fLL,fUL] = falltime(allsums);    flight_ends = round(fLT+tracking_Fs/2);       %... +Fs is a sistematic correction, useful
num_flights = size(R,2);
ref = ones(size(flight_starts));
avg_flight_time = mean((flight_ends-flight_starts)./tracking_Fs);

%% Plot results

% Plot sampling interval counts
if analyze_Ca
    figure();
    hist(img_sampling); set(gca, 'YScale', 'log'); ylabel('counts');    xlabel('interframe(s)');
end

% Plot 2D flight trajectories
figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
plot(x_mean(1,:),x_mean(2,:),'.');
hold on;        rectangle('Position',[xL yB xR-xL yF-yB]);
scatter([F3(1) F4(1)],[F3(2) F4(2)],'filled');  hold off;
xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);

subplot(1,2,2);
plot(x_spl(1,:),x_spl(2,:)); hold on; plot(x_mean(1,:),x_mean(2,:),'.','MarkerSize',1);
rectangle('Position',[xL yB xR-xL yF-yB]);  hold off;
xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);

% Plot session timeline
figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
ax1 = subplot(3,1,1);   plot(t,x_mean(1,:),'.');  hold on;
plot(new_t,x_spl(1,:),'--');   refline(0,F3(1));    hold off;
legend('cluster/mean','spl');   ylabel('x (m)');
ax2 = subplot(3,1,2);   plot(new_t,v_abs,'.');
hold on;    stem(new_t(flight_starts),ref);    stem(new_t(flight_ends),ref);  hold off;
ylabel('v (m/s)');
ax3 = subplot(3,1,3);   plot(t,rew_signal);
ylabel('Rewards');
linkaxes([ax1,ax2,ax3],'x');    xlabel('Samples');

%% Clustering flights
%Cut out flights, downsample to ds_clus positions per flight
all_flights = NaN(3,max(flight_ends-flight_starts),num_flights);    %3D matrix with all flights
all_flights_ds = NaN(3,ds_clus,num_flights);                        %3D matrix with all flights(downsampled)

for nf = 1 : size(all_flights,3)
    trajectory = x_spl(:,flight_starts(nf):flight_ends(nf));
    velocity = v_abs(:,flight_starts(nf):flight_ends(nf));
    all_flights(:,1:(flight_ends(nf)-flight_starts(nf))+1,nf) = trajectory;
    all_flights_vel(1,1:(flight_ends(nf)-flight_starts(nf)+1),nf) = velocity;
    all_flights_ds(:,:,nf) = interp1(linspace(1,3,size(trajectory,2)),trajectory',linspace(1,3,ds_clus),'spline')';
    
    %     %Uncomment if you want to see how the downsampled flights look like
    %     plot3(all_flights(1,:,nf),all_flights(2,:,nf),all_flights(3,:,nf),'Color','b');
    %     hold on;
    %     plot3(all_flights_ds(1,:,nf),all_flights_ds(2,:,nf),all_flights_ds(3,:,nf),'Color','r');
    %     hold off;
    %     w = waitforbuttonpress;
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
    figure();
    Y = pdist(X,'euclidean');   Z = linkage(Y);
    hLines = dendrogram(Z,0);  hold on;    refline(0,dist);    hold off;
    idx = cluster(Z,'Cutoff',dist,'Criterion','distance');
    title([num2str(length(unique(idx))) ' clusters']);
    ylim([0 10]);
end

% Create structure with flight start stop frames, id of the trajectory
clear flight;
flight.strt_frame = ceil(flight_starts)';
flight.stop_frame = ceil(flight_ends)';
flight.pos = all_flights;
flight.vel = all_flights_vel;
flight.id = idx;
flight.Fs = tracking_Fs;

% Sort structure according to cluster id
clear flight_sorted;
[flight_sorted.id,I] = sort(flight.id);
flight_sorted.strt_frame = flight.strt_frame(I);
flight_sorted.stop_frame = flight.stop_frame(I);
flight_sorted.pos = flight.pos(:,:,I);
flight_sorted.vel = flight.vel(:,:,I);
flight_sorted.Fs = flight.Fs;
flight_sorted.N = size(flight_sorted.id,1);

% Assign isolated clusters to cluster #flights+1
[Ns,b] = histc(flight_sorted.id,unique(flight_sorted.id));
flight_sorted.id(Ns(b)<N_min) = size(all_flights,3)+1;             %flight_sorted.id(Ns(b)==1) = size(all_flights,3)+1;
id_surv_clusters = unique(flight_sorted.id);
n_surv_clusters = size(id_surv_clusters,1);

% Create final structure flight.clus after re-assignment
clear flight_clus;
flight_clus.id = flight_sorted.id;
flight_clus.strt_frame = flight_sorted.strt_frame;
flight_clus.stop_frame = flight_sorted.stop_frame;
flight_clus.pos = flight_sorted.pos;
flight_clus.vel = flight_sorted.vel;
flight_clus.Fs = flight_sorted.Fs;
flight_clus.N = flight_sorted.N;
for jj=1:n_surv_clusters;
    flight_clus.id(flight_sorted.id == id_surv_clusters(jj)) = jj;
end
id_surv_clusters = unique(flight_clus.id);

%Re-assign id for convenience, if necessary
if reassign
    new_ord = [];
    [~,new_ord] = sort(histc(flight_clus.id,id_surv_clusters(1:end-1)),'descend');
    new_ord = [new_ord; id_surv_clusters(end)];
    new_ord = circshift(new_ord,1);
    reassign_matrix =(flight_clus.id == new_ord');
    for jj=1:n_surv_clusters;
        flight_clus.id(reassign_matrix(:,jj)) = jj;
    end  
end 

%Calculate trajectory length, duration in s and interflight (take-off to take-off)
for ii = 1:flight_sorted.N
    flight_clus.length(ii)= arclength(flight_clus.pos(1,~isnan(flight_clus.pos(1,:,ii)),ii),flight_clus.pos(2,~isnan(flight_clus.pos(2,:,ii)),ii),flight_clus.pos(3,~isnan(flight_clus.pos(3,:,ii)),ii),'s');
    flight_clus.dur(ii) = (flight_clus.stop_frame(ii)-flight_clus.strt_frame(ii))./flight_clus.Fs;
end
flight_clus.ifd = diff(flight_clus.strt_frame)';

%Visualize
figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
col = hsv(n_surv_clusters);
for jj=1:n_surv_clusters;
    id = find(flight_clus.id==jj);
    
    subplot(3,n_surv_clusters,jj);
    avg_take_off = [];
    for ii=1:size(id,1);
        hold on;
        title(['Cluster' num2str(jj) '  (' num2str(size(id,1)) ' flights)'])
        plot3(flight_clus.pos(1,:,id(ii)),flight_clus.pos(2,:,id(ii)),flight_clus.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col(jj,:));
        avg_take_off = [avg_take_off flight_clus.pos(:,1,id(ii))];
        hold on;
    end
    take_off = mean(avg_take_off,2);
    textscatter3(take_off(1),take_off(2),take_off(3),"Take-off");
    
    plot3(x_spl(1,:),x_spl(2,:),x_spl(3,:),':','Color',[0.7 0.7 0.7],'MarkerSize',0.001);
    xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
    xlabel('x');    ylabel('y');    zlabel('z');    view(2);
    hold off;
    
    subplot(3,n_surv_clusters,n_surv_clusters+jj);
    avg_take_off = [];
    for ii=1:size(id,1);
        hold on;
        title(['Cluster' num2str(jj) '  (' num2str(size(id,1)) ' flights)'])
        plot3(flight_clus.pos(1,:,id(ii)),flight_clus.pos(2,:,id(ii)),flight_clus.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col(jj,:));
        avg_take_off = [avg_take_off flight_clus.pos(:,1,id(ii))];
        hold on;
    end
    take_off = mean(avg_take_off,2);
    textscatter3(take_off(1),take_off(2),take_off(3),"Take-off");
    
    plot3(x_spl(1,:),x_spl(2,:),x_spl(3,:),':','Color',[0.7 0.7 0.7],'MarkerSize',0.001);
    xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
    xlabel('x');    ylabel('y');    zlabel('z');    view(0,0);
    hold off;
    
    subplot(3,n_surv_clusters,n_surv_clusters*2+jj);
    histogram(flight_clus.dur(id));
    xlim([0 15]);   xlabel('Duration(s)');  ylabel('Counts');
end