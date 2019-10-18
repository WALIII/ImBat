function [CutCells, rhoAll, Template] = ImBat_ClusterDemon(ROI_Data,out,day,clust);

% 'out' variable is from [out] =  ImBat_SegTrajectories
% WAL3
% 08/08/2019

% Variables:
ROI_ON = 100;
ROI_OFF = 400;

% Pad data:
ForePad = 10;
AftPad = 300;


% Definitions:
A = ROI_Data{1,day}.Alignment.out.flights; % flight Data
At = ROI_Data{1,day}.Alignment.out.Location_time;
Ct = ROI_Data{1,day}.Alignment.out.video_times; % flight times

% % Cut out flights
b = isnan(A(:,1));
b = double(-b+1);
[w1 w2] = findpeaks(b);

idX = out.ClusterIndex{clust};
for ii = 1:size(idX,2)
    ClustFlight(:,:,ii) = A(out.flight_starts_indx(idX(ii))-ForePad:out.flight_starts_indx(idX(ii))+AftPad,:);
    InFlight(:,ii) = b(out.flight_starts_indx(idX(ii))-ForePad:out.flight_starts_indx(idX(ii))+AftPad);
    FlightTimes(:,ii) = At(out.flight_starts_indx(idX(ii))-ForePad:out.flight_starts_indx(idX(ii))+AftPad);
end

figure();
hold on;
for i = 1:size(ClustFlight,3)
    plot(squeeze(ClustFlight(:,1,i)),'r')
    plot(squeeze(ClustFlight(:,2,i)),'g')
    plot(squeeze(ClustFlight(:,3,i)),'b')
end


%% Cut out Ca Imaging data for each flight:

% find closest start/end times
for i = 1: size(FlightTimes,2);
    [minValue_1(:,i),closestIndex_1(:,i)] = min(abs(FlightTimes(1,i)-Ct));
    [minValue_2(:,i),closestIndex_2(:,i)] = min(abs(FlightTimes(size(FlightTimes,1),i)-Ct)); % for posterity
end
% typical length:
RoundFrames = round((FlightTimes(size(FlightTimes,1)) - FlightTimes(1,1) )*30);


% plot individual cells
roidat = ROI_Data{1,day}.ROIs.results.C;
roidat2 = (ROI_Data{1,day}.ROIs.results.C_raw);

% Frame offsets ( a new pad...);
for i = 1:size(FlightTimes,2);
    try   
        CutCells(:,:,i) = roidat(:,closestIndex_1(:,i)-ROI_ON :closestIndex_1(:,i)+RoundFrames+ROI_OFF);
        CutCells2(:,:,i) = roidat2(:,closestIndex_1(:,i)-ROI_ON :closestIndex_1(:,i)+RoundFrames+ROI_OFF);
        Ydata(:,:,:,i) = Y(:,:,closestIndex_1(:,i)-ROI_ON :closestIndex_1(:,i)+RoundFrames+ROI_OFF);
    catch
        disp('Flight is too close to beginning, or end of flight for current pad setting');
    end
end



% index of the most likely 'place cells'
mC = squeeze(mean(CutCells(:,:,:),3));
mCA = (squeeze(mean(CutCells(:,:,1:2:size(CutCells,3)),3)));
mCB = (squeeze(mean(CutCells(:,:,2:2:size(CutCells,3)),3)));
figure();
for i = 1:size(mCA,1)
    mCD(i) = corr(zscore(mCA(i,30:400)'),zscore(mCB(i,30:400)'),'type','Spearman');
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
for ii = 1:size(mC,1);
    i = mc42(ii);
    Dplot = squeeze(CutCells(i,:,:))';
    imagesc(zscore(Dplot,[],2));
    pause(0.01);
end



clear temp data
for i = 1:size(CutCells,3);
    temp(i,:,:) = CutCells(:,50:500,i)';
end
data.undirected = temp(1:2:size(CutCells,3),:,:);
data.directed = temp(2:2:size(CutCells,3),:,:);
figure();
[indX99,B,C, output] = CaBMI_schnitz(data);

for ii = 1:size(output.Index,1)
    ix = output.Index(ii);
    [rho,pval] = corr(B(ix,:)',C(ix,:)');
    rhoAll(:,ii) = rho;
end

Template = mean(CutCells(:,50:250,:),3);

%Template = output.FullSort(output.Index,40:200);

disp([num2str(size(CutCells,3)), ' Flights']);


