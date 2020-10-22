function ImBat_Align_FC(CombinedROI,flightPaths);
% Align Calcium and FLight data across the interval



% First, put all aligned time-series data in the same place, it should be
% on the same clock
AlignedData.C = 
AlignedData.C_r = 
AlignedData.S = 
AlignedData.Flights = flightPaths.AllFlights



% Input paramaters:
ForePad = 
AftPad = 


clust = 2;

A = flightPaths.AllFlights;
At = flightPaths.AllFlightsMasterTime;
Velocity = flightPaths.batSpeed;
Ct

for i = clust; % for this clustered Trajectory ( use one to start)
    idX = out.ClusterIndex{i};
    for ii = 1:size(idX,2)
        try
            ClustFlight(:,:,ii) = A(flightPaths.flight_starts_indx(idX(ii))-ForePad:flightPaths.flight_starts_indx(idX(ii))+AftPad,:);
            FlightTimes(:,ii) = At(flightPaths.flight_starts_indx(idX(ii))-ForePad:flightPaths.flight_starts_indx(idX(ii))+AftPad);
%             CutFlights(:,ii) = Velocity(flightPaths.flight_starts_indx(idX(ii))-(ROI_ON/fs)*tfs:flightPaths.flight_starts_indx(idX(ii))+(ROI_OFF/fs)*tfs );
        catch
            disp('flight too close to end')
        end
    end
end



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







% 
% try
%         ClustFlight{i}(:,:,ii) = Flights(flightPaths.flight_starts_indx(idX{i}(ii))-ForePad:out.flight_starts_indx(idX{i}(ii))+AftPad,:);
%         InFlight{i}(:,ii) = b(out.flight_starts_indx(idX{i}(ii))-ForePad:out.flight_starts_indx(idX{i}(ii))+AftPad);
%         FlightTimes{i}(:,ii) = Location_time(out.flight_starts_indx(idX{i}(ii))-ForePad:out.flight_starts_indx(idX{i}(ii))+AftPad);
%         CutFlights{i}(:,ii) = Velocity(out.flight_starts_indx(idX{i}(ii))-(ROI_ON/fs)*tfs:out.flight_starts_indx(idX{i}(ii))+(ROI_OFF/fs)*tfs );
%         catch
%             disp('flight too close to end');
%         end
%         end
%     % Get Calcium Data:
%     %% Cut out data for each flight:
% 
%     % find closest start/end times
%     for iii = 1: size(FlightTimes{i},2);
%         [minValue_1{i}(:,iii),closestIndex_1{i}(:,iii)] = min(abs(FlightTimes{i}(1,iii)-Video_times));
%         [minValue_2{i}(:,iii),closestIndex_2{i}(:,iii)] = min(abs(FlightTimes{i}(size(FlightTimes{i},1),iii)-Video_times)); % for posterity
%     end
%     % typical length:
%     RoundFrames = round((FlightTimes{i}(size(FlightTimes{i},1)) - FlightTimes{i}(1,1) )*fs);
%     %RoundFrames = 0;
% 
%     for iii = 1:size(FlightTimes{i},2);
%         try
% 
%             CutCells{i}(:,:,iii) = roidat(:,closestIndex_1{i}(:,iii):closestIndex_1{i}(:,iii)+RoundFrames);
%             CutCells2{i}(:,:,iii) = roidat2(:,closestIndex_1{i}(:,iii) :closestIndex_1{i}(:,iii)+RoundFrames);
%             if exist('Y')
%                 Ydata{i}(:,:,:,i) = Y(:,:,closestIndex_1{i}(:,i)-ROI_ON :closestIndex_1{i}(:,i)+RoundFrames+ROI_OFF);
%             end
%         catch
%             disp('Flight is too close to beginning, or end of flight for current pad setting');
%         end
%     end
% 