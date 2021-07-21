flight_type_ifd = 2;

true_times = flightPaths34.AllFlightsMasterTime(flightPaths34.flight_starts_idx(:));
[tts,tts_idx]  = sort(true_times(:));
%[r_tts,r_tts_idx] = sort(reward_times(:));
dayx_flight_timeline = flightPaths34.flight_starts_idx(tts_idx)';

temp_F = find(flightPaths34.id==flight_type_ifd);
temp_F = temp_F(1:end-1);
temp_ifd = flightPaths34.ifd(temp_F);
IFD_list = temp_ifd(temp_ifd>0);

figure(); histogram(IFD_list);

% Remove outliers
IFD_noOutliers = rmoutliers(IFD_list);

figure(); histogram(IFD_noOutliers);
outlier_edges = [min(IFD_noOutliers),max(IFD_noOutliers)];

IFD_Megalist = {};

%% Bin by day to see if the IFD for given flights gets tighter
for j=1:max(fd)
    day1_indexs = find(fd==j);
    day1_flight_timeline = dayx_flight_timeline(find(fd==j));
    day_1_flights = c_s_34(find(fd==j));
    day_1_ifd = diff(day1_flight_timeline)';
    
    temp_F = find(day_1_flights==flight_type_ifd);
    temp_F = temp_F(1:end-1);
    if isempty(temp_F)
        disp(strcat("No flight type ",num2str(flight_type_ifd)," on day ",num2str(j)));
        continue
    end
    temp_ifd = day_1_ifd(temp_F);
    IFD_list = temp_ifd(temp_ifd>0);
    
    IFD_Megalist{j} = IFD_list;
    figure(); hold on; histogram(IFD_list); title(strcat("Day ",num2str(j)," inter-flight interval of all ",num2str(flight_type_ifd), " flight types. (Total ",num2str(size(temp_F,1)),")"));
    
    IFD_means(j) = mean(IFD_list);
    IFD_vars(j) = var(IFD_list);
end

figure(); hold on; plot(nonzeros(IFD_means)); title(strcat("Mean Interflight Interval. Flight type ",num2str(flight_type_ifd)));
figure(); hold on; plot(nonzeros(IFD_vars)); title(strcat("Variance of Interflight Interval. Flight type ",num2str(flight_type_ifd)));

for i=1:size(IFD_Megalist,2)
    IFD_mean(i) = mean(IFD_Megalist{i});
end

figure(); bar(IFD_mean); ylabel('Interflight Interval (ms)');

%% Identify places where there are more than 4 in a row and see if over time the IFI gets smaller
test_seq = [5 2 2 2 2]';
test_seq = [2 2 2 2 2]';
test_seq = [5 2]';
test_seq = [2 3 3]';
test_seq = [1 1 1]';
test_seq = [2 2 2]';
test_seq = [10 1]';
test_seq = [1 4]';
test_seq = [4 4 4]';
test_seq = [1 6 1]';
test_seq = [1 1 1]';
test_seq = [1 5 5]';

RDS = [];

lim = size(test_seq,1)-1;
day_1_ifd_seqs = [];
for j=1:max(fd)
    day1_indexs = find(fd==j);
    day1_flight_timeline = dayx_flight_timeline(find(fd==j));
    day_1_flights = c_s_34(find(fd==j));
    for r=1:size(day_1_flights,1)-lim
        sseq = day_1_flights(r:r+lim);
        if sseq == test_seq
            day_1_ifd_seqs = [day_1_ifd_seqs;diff(day1_flight_timeline(r:r+lim))'];
        end
    end 
    ds_day_1_ifd_seqs = day_1_ifd_seqs(1:lim:end,:);
    rds_day_1_ifd_seqs = reshape(ds_day_1_ifd_seqs,1,size(ds_day_1_ifd_seqs,1)*size(ds_day_1_ifd_seqs,2));
    %rds_day_1_ifd_seqs = rds_day_1_ifd_seqs(rds_day_1_ifd_seqs<10000);

    RDS = [RDS, rds_day_1_ifd_seqs];
%     figure(); hold on;
%     plot(rds_day_1_ifd_seqs);
%     yline(mean(rds_day_1_ifd_seqs));
%     disp(size(rds_day_1_ifd_seqs,2));
%     disp(var(rds_day_1_ifd_seqs));

    x = [1:size(rds_day_1_ifd_seqs,2)];
    y = rds_day_1_ifd_seqs;
    x1 = [1:size(rds_day_1_ifd_seqs,2)];
    poly = polyfit([1:size(rds_day_1_ifd_seqs,2)],rds_day_1_ifd_seqs,1);
    y1 = polyval(poly,x1);
%     figure();
%     plot(x,y,'o')
%     hold on;
%     plot(x1,y1)
%     hold off;
end

x = [1:size(RDS,2)];
y = RDS;
x1 = [1:size(RDS,2)];
poly = polyfit([1:size(RDS,2)],RDS,1);
y1 = polyval(poly,x1);
figure();
plot(x,y,'o')
hold on;
plot(x1,y1)
hold off;

figure(); hold on;
plot(RDS);
yline(mean(RDS));
%disp(size(RDS,2));
%disp(var(RDS));




% ds_day_1_ifd_seqs = day_1_ifd_seqs(1:lim:end,:);
% rds_day_1_ifd_seqs = reshape(ds_day_1_ifd_seqs,1,size(ds_day_1_ifd_seqs,1)*size(ds_day_1_ifd_seqs,2));
% rds_day_1_ifd_seqs = rds_day_1_ifd_seqs(rds_day_1_ifd_seqs<10000);
% 
% figure(); hold on;
% plot(rds_day_1_ifd_seqs);
% yline(mean(rds_day_1_ifd_seqs));
% disp(size(rds_day_1_ifd_seqs,2));
% disp(var(rds_day_1_ifd_seqs));
% 
% x = [1:size(rds_day_1_ifd_seqs,2)];
% y = rds_day_1_ifd_seqs;
% x1 = [1:size(rds_day_1_ifd_seqs,2)];
% poly = polyfit([1:size(rds_day_1_ifd_seqs,2)],rds_day_1_ifd_seqs,1);
% y1 = polyval(poly,x1);
% figure();
% plot(x,y,'o')
% hold on;
% plot(x1,y1)
% hold off;