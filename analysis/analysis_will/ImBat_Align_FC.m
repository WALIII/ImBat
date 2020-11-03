function [FlightAlignedROI] = ImBat_Align_FC(CombinedROI,flightPaths,clust);
% Align Calcium and FLight data across the interval


fs = 30;
tfs = 120;

% User paramaters:
ForePad = 500;
AftPad = 2500;
ROI_ON = (ForePad/tfs)*fs;
ROI_OFF = (AftPad/tfs)*fs;




A = flightPaths.AllFlights;
At = flightPaths.AllFlightsMasterTime;
Velocity = flightPaths.batSpeed;
Ct = CombinedROI.timestamps';

for i = clust; % for this clustered Trajectory ( use one to start)
    idX = flightPaths.clusterIndex{i}';
    for ii = 1:size(idX,2)
        try
            ClustFlight(:,:,ii) = A(flightPaths.flight_starts_idx(idX(ii))-ForePad:flightPaths.flight_starts_idx(idX(ii))+AftPad,:);
            FlightTimes(:,ii) = At(flightPaths.flight_starts_idx(idX(ii)));
            FlightLength(:,ii) = flightPaths.flight_ends_idx(idX(ii))-flightPaths.flight_starts_idx(idX(ii));
            %             CutFlights(:,ii) = Velocity(flightPaths.flight_starts_indx(idX(ii))-(ROI_ON/fs)*tfs:flightPaths.flight_starts_indx(idX(ii))+(ROI_OFF/fs)*tfs );
        catch
            disp('flight too close to end')
        end
    end
end



figure();
hold on;
for i = 1:size(ClustFlight,3)
    plot3(squeeze(ClustFlight(ForePad-40:end-AftPad/2,1,i)),squeeze(ClustFlight(ForePad-40:end-AftPad/2,2,i)),squeeze(ClustFlight(ForePad-40:end-AftPad/2,3,i)),'r')
end
grid on;

figure();
hold on;
for i = 1:size(ClustFlight,3)
    plot(((1:size(ClustFlight,1))/120)-ForePad/tfs,squeeze(ClustFlight(:,1,i)),'r')
    plot(((1:size(ClustFlight,1))/120)-ForePad/tfs,squeeze(ClustFlight(:,2,i)),'g')
    plot(((1:size(ClustFlight,1))/120)-ForePad/tfs,squeeze(ClustFlight(:,3,i)),'b')
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
roidat = CombinedROI.C;
roidat2 = CombinedROI.C_raw;
roidat3 = full(CombinedROI.S);

% smooth data
disp('Smoothing data...');
for i = 1:size(roidat2,1)
    roidat(i,:) = smooth(roidat(i,:),30);
    roidat2(i,:) = smooth(roidat2(i,:),15);
    roidat3(i,:) = smooth(roidat3(i,:),30);
end


for i = 1:size(FlightTimes,2);
    try
        
        CutCells(:,:,i) = roidat(:,closestIndex_1(:,i)-ROI_ON :closestIndex_1(:,i)+ROI_OFF);
        CutCells2(:,:,i) = roidat2(:,closestIndex_1(:,i)-ROI_ON :closestIndex_1(:,i)+ROI_OFF);
        CutCells3(:,:,i) = roidat3(:,closestIndex_1(:,i)-ROI_ON :closestIndex_1(:,i)+ROI_OFF);
        CutCells_date(i) = CombinedROI.day_vector(closestIndex_1(:,i));
        if exist('Y')
            Ydata(:,:,:,i) = Y(:,:,closestIndex_1(:,i)-ROI_ON :closestIndex_1(:,i)+ROI_OFF);
        end
    catch
        disp('Flight is too close to beginning, or end of flight for current pad setting');
    end
end


% figure();
% for i = 1:size(CutCells,1);
%     adata = zscore(squeeze(CutCells2(i,:,:)),[],1)';
%     imagesc(adata);
% % Get axis handle
% ax = gca;
% % Set where ticks will be
% ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
% % Set TickLabels;
% ylabel('filghts')
% xlabel('time from takeoff');
% ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};
%     pause();
% end

col = jet(90);

figure();
CutCells_nan = CutCells;
CutCells_nan(CutCells_nan==0) = NaN;
hold on;
for i = 1:size(CutCells,1);;
adata = zscore(squeeze(CutCells_nan(i,:,:)),[],1)';
mean_dat(:,i) = nanmean(adata);
adata = zscore(squeeze(CutCells(i,:,:)),[],1)';
adata = adata+2*i;

L = size(adata,2);
se = std(adata)/2;%/10;%sqrt(length(adata));
mn = nanmean(adata);
h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col(i,:)); alpha(0.5);
plot(mn,'Color',col(i,:));
end

mean_dat(isnan(mean_dat))=0;
% Sort data
[aa ab] = max(mean_dat(100:400,:));

thresh2use = 0.75;
indx2sort = find(aa>thresh2use);% find peaks greater than 1 and sort these ones:
index_rest = find(aa<=thresh2use);
first_idx = cat(2,indx2sort,index_rest);
covert1 = mean_dat(:,first_idx);
ab = ab(first_idx);

figure(); imagesc(covert1');

[a2 b2] = sort(ab(1:size(indx2sort,2)),'ascend');

index_rest2 = max(b2)+1:size(index_rest,2)+max(b2);
final_idx = cat(2,b2,index_rest2);

covert2 = covert1(:,final_idx);

figure(); imagesc(covert2');
% now create the index:

% create index to export:
clear IDX;
IDX = 1:size(covert2,2);
IDX = IDX(first_idx);
IDX = IDX(final_idx);

figure(); imagesc(mean_dat(:,IDX)');

FlightAlignedROI.C = CutCells;
FlightAlignedROI.C_raw = CutCells2;
FlightAlignedROI.S = CutCells3;
FlightAlignedROI.clust_number = clust;
FlightAlignedROI.IDX = IDX;
FlightAlignedROI.ClustFlight_withPads = ClustFlight;
FlightAlignedROI.ClustFlight = ClustFlight(ForePad-40:end-AftPad/2,:,:);
FlightAlignedROI.FlightLength = FlightLength;

FlightAlignedROI.ForePad = ForePad;
FlightAlignedROI.AftPad = AftPad;
FlightAlignedROI.ROI_ON = ROI_ON;
FlightAlignedROI.ROI_OFF = ROI_OFF;
FlightAlignedROI.CutCells_date =  CutCells_date;



