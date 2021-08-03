function [Q1 Q4] = ImBat_BehavStabilityVsTime(out_markov,flightPaths,flight2use)
% look at the stability of behavior over time



idx2use = find(out_markov.FlightIDVector == flight2use);
% build flight paths:
for i = 1:size(flightPaths.flight_starts_idx,2) 
try
    Flights2Use(:,:,i) =  flightPaths.AllFlights(flightPaths.flight_starts_idx(i)-120:flightPaths.flight_starts_idx(i)+(10*120),:);
catch
    Flights2Use(:,:,i) =  flightPaths.AllFlights(flightPaths.flight_starts_idx(i-1)-120:flightPaths.flight_starts_idx(i-1)+(10*120),:);
end
end

Flights2Use = Flights2Use(:,:,out_markov.out_sort); % get the right flights in the right order
Flights2Use = Flights2Use(:,:,idx2use);


% check flights:
% figure(); 
% tidx = (1:length(squeeze(Flights2Use(:,1,1))))/120;
% for i = 1:3;
% subplot(1,3,i);
% hold on;
% plot(tidx,squeeze(Flights2Use(:,i,:)),'b');
% end

% now, take the corr:
Mflight = nanmedian(Flights2Use,3); % mean flight

% get XYZ


% get offset:
endOffsets = (flightPaths.flight_ends_idx-flightPaths.flight_starts_idx);
endOffsets = endOffsets(out_markov.out_sort);
endOffset = round(mean(endOffsets(idx2use)))+240;


for i = 1:size(Flights2Use,3)
tXYZ(i,1) = corr(Mflight(1:endOffset,1),squeeze(Flights2Use(1:endOffset,1,i)));
tXYZ(i,2) = corr(Mflight(1:endOffset,2),squeeze(Flights2Use(1:endOffset,2,i)));
tXYZ(i,3) = corr(Mflight(1:endOffset,3),squeeze(Flights2Use(1:endOffset,3,i)));
end

% figure(); plot((tXYZ')');

% figure(); 
hold on;
plot(mean(tXYZ')','.');
% plot(smooth(mean(tXYZ'),50),'--','LineWidth',3);

GG = mean(tXYZ')';

% first and last quarter:
Q1 = mean(GG(1:round(size(GG,1)/4)));
Q4 = mean(GG((round(size(GG,1)/4))*2:round(size(GG,1)/4)*3-4));

