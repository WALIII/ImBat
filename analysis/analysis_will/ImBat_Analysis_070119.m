function [CutCells, Ydata, ClustFlight,CutCells2, CutFlights,Velocity] = ImBat_Analysis_070119(ROI_Data,varargin);

%% align ROI data to the cut out flight data...

% Which day
day = 1;
clust = 1;
LoadY =0;
fs = 30; % video fs
tfs = 120; % tracking fs
% Frame offsets ( a new pad...);
% Frame offsets ( a new pad...);
ROI_ON = 100;
ROI_OFF = 100;
% Pad data:
ForePad = (ROI_ON/fs)*tfs;
AftPad = (ROI_OFF/fs)*tfs;
% Manual inputs
vin=varargin;
for i=1:length(vin)
    if isequal(vin{i},'clust') % manually inputing a sort order
        clust=vin{i+1};
    elseif isequal(vin{i},'day')
        day=vin{i+1};
        fs = ROI_Data{day}.ROIs.results.Fs; % video fs
        
    elseif isequal(vin{i},'HL')
        HL = vin{i+1};
    elseif isequal(vin{i},'Y')
        Y = vin{i+1};
    end
    
end

if LoadY ==1;
    if exist('Y') ==0; % load in Y from local directory
        disp( 'Y matrix is being loaded from local directory...');
        load([ROI_Data{day}.date,'/',ROI_Data{day}.folder,'/Motion_corrected_Data_DS.mat'])
    end
end

disp(['Running day ', ROI_Data{day}.date]);
% Plot flights
figure();
a(1) = subplot(2,1,1);
hold on
plot(ROI_Data{1,day}.Alignment.out.Location_time,(ROI_Data{1,day}.Alignment.out.flights)./500);
len = size(ROI_Data{1,day}.ROIs.results.S,2)-5;;
plot(ROI_Data{1,day}.Alignment.out.video_times2(1:len),zscore(smooth(mean(full(ROI_Data{1,day}.ROIs.results.S(1:10,1:len)),1),100)));
a(2) = subplot(2,1,2);
hold on;
% C = mean(full(ROI_Data{1,day}.ROIs.results.C(1:50,:)));
Ct = ROI_Data{1,day}.Alignment.out.video_times2;
num2plot = size(ROI_Data{1,1}.ROIs.results.C,1);

if num2plot>50; % only plot 50
    num2plot = 50;
end

for i = 1: 10;
    C1 = zscore(ROI_Data{1,day}.ROIs.results.C(i,:))+i*2;
    plot(Ct, C1);
end

linkaxes(a, 'x');


%% cluster flights
A = ROI_Data{1,day}.Alignment.out.flights; % flight Data
At = ROI_Data{1,day}.Alignment.out.Location_time;
Ct = ROI_Data{1,day}.Alignment.out.video_times2; % flight times

% % Cut out flights
b = isnan(A(:,1));
b = double(-b+1);
[w1 w2] = findpeaks(b);

[out] =  ImBat_SegTrajectories(ROI_Data{1,day}.Alignment.out.flights,ROI_Data{1,day}.Alignment.out.Location_time);


% get velocity:
Velocity = (abs(diff(ROI_Data{1,day}.Alignment.out.Location2(:,1)))+abs(diff(ROI_Data{1,day}.Alignment.out.Location2(:,2)))+abs(diff(ROI_Data{1,day}.Alignment.out.Location2(:,3))))/3;
Velocity = smooth(Velocity,100);

for i = clust; % for this clustered Trajectory ( use one to start)
    idX = out.ClusterIndex{i};
    for ii = 1:size(idX,2)
        try
            ClustFlight(:,:,ii) = A(out.flight_starts_indx(idX(ii))-ForePad:out.flight_starts_indx(idX(ii))+AftPad,:);
            InFlight(:,ii) = b(out.flight_starts_indx(idX(ii))-ForePad:out.flight_starts_indx(idX(ii))+AftPad);
            FlightTimes(:,ii) = At(out.flight_starts_indx(idX(ii))-ForePad:out.flight_starts_indx(idX(ii))+AftPad);
            CutFlights(:,ii) = Velocity(out.flight_starts_indx(idX(ii))-(ROI_ON/fs)*tfs:out.flight_starts_indx(idX(ii))+(ROI_OFF/fs)*tfs );
        catch
            disp('flight too close to end')
        end
    end
end
figure(); imagesc(InFlight')

figure();
hold on;
for i = 1:size(ClustFlight,3)
    plot3(squeeze(ClustFlight(:,1,i)),squeeze(ClustFlight(:,2,i)),squeeze(ClustFlight(:,3,i)),'r')
end

figure();
hold on;
for i = 1:size(ClustFlight,3)
    plot(squeeze(ClustFlight(:,1,i)),'r')
    plot(squeeze(ClustFlight(:,2,i)),'g')
    plot(squeeze(ClustFlight(:,3,i)),'b')
end



%% Cut out data for each flight:

% find closest start/end times
for i = 1: size(FlightTimes,2);
    [minValue_1(:,i),closestIndex_1(:,i)] = min(abs(FlightTimes(1,i)-Ct));
    [minValue_2(:,i),closestIndex_2(:,i)] = min(abs(FlightTimes(size(FlightTimes,1),i)-Ct)); % for posterity
end
% typical length:
RoundFrames = round((FlightTimes(size(FlightTimes,1)) - FlightTimes(1,1) )*fs);

% plot individual cells
roidat = ROI_Data{1,day}.ROIs.results.C;
roidat2 = (ROI_Data{1,day}.ROIs.results.C_raw);

% smooth data
disp('Smoothing data...');
for i = 1:size(roidat2,1)
    roidat2(i,:) = smooth(roidat2(i,:),15);
end


for i = 1:size(FlightTimes,2);
    try
        
        CutCells(:,:,i) = roidat(:,closestIndex_1(:,i)-ROI_ON :closestIndex_1(:,i)+RoundFrames+ROI_OFF);
        CutCells2(:,:,i) = roidat2(:,closestIndex_1(:,i)-ROI_ON :closestIndex_1(:,i)+RoundFrames+ROI_OFF);
        if exist('Y')
            Ydata(:,:,:,i) = Y(:,:,closestIndex_1(:,i)-ROI_ON :closestIndex_1(:,i)+RoundFrames+ROI_OFF);
        end
    catch
        disp('Flight is too close to beginning, or end of flight for current pad setting');
    end
end



% index of the most likely 'place cells'

% mC = squeeze(mean(CutCells(:,:,:),3));
% mCA = (squeeze(mean(CutCells(:,:,1:2:size(CutCells,3)),3)));
% mCB = (squeeze(mean(CutCells(:,:,2:2:size(CutCells,3)),3)));
% figure();
% for i = 1:size(mCA,1)
% mCD(i) = corr(zscore(mCA(i,30:400)'),zscore(mCB(i,30:400)'),'type','Spearman');
% end

for i = 1: size(CutCells,1)
    a = corr(squeeze(CutCells(i,:,:)));
    mCD(:,i) = median(mean((a)));
end
% hold on;
% plot(mCA(i,:));
% plot(mCB(i,:));
% hold off
% title([ 'Pearson Correlation of ', num2str(mCD(i))]);
%
% pause();
% clf
% end


% mC2 = (max(mC(:,40:300),[],2)-min(mC(:,40:300),[],2))./(max(mC(:,400:500),[],2)-min(mC(:,400:500),[],2));

% [mc31, mc32] = sort(mC2,'descend');
[mc41, mc42] = sort(mCD,'descend');


figure();
histogram(mCD);
figure();
for ii = 1:size(mc42,2);
    i = mc42(ii);
    Dplot = squeeze(CutCells(i,:,:))';
    imagesc(zscore(Dplot,[],2));
    title(['Cell number ' num2str(i)]);
    pause(0.01);
end



clear temp data
for i = 1:size(CutCells,3);
    temp(i,:,:) = CutCells(:,:,i)';
end
data.undirected = temp(1:2:size(CutCells,3),:,:);
data.directed = temp(2:2:size(CutCells,3),:,:);
figure();
[indX99,B,C, output] = CaBMI_schnitz(data);
disp([num2str(size(CutCells,3)), ' Flights']);

% Now, largely cluster:
roidat_2 = zscore(roidat(indX99,:),[],2);

numFlights = 7;
idx = kmeans(roidat_2,numFlights);
[~,idx2] = sort(idx);
figure(); imagesc(roidat_2(idx2,:),[0 10]);
colormap(hot);


if exist('Y')
else
    Ydata = 0;
end
end
%
%
% figure(); imagesc(roidat2(indX99,:),[0 10]);
%
%
%
% figure();
% col = lines(120);
% hold on;
% for i = 61:120;
%     for ii = 1:size(CutCells2,3);
% Dplot = zscore(squeeze(CutCells2(mc32(i),:,ii)),[],2)+i*8;;
% plot(smooth(Dplot,10),'Color',col(i,:));
%     end
% end
%
%

