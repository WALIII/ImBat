function [Q1 Q4 out] = ImBat_BehavStabilityVsTime(out_markov,flightPaths,flight2use)
% look at the stability of behavior over time



idx2use = find(out_markov.FlightIDVector == flight2use);
% build flight paths:
for i = 1:size(flightPaths.flight_starts_idx,2)
    try
        Flights2Use(:,:,i) =  flightPaths.AllFlights(flightPaths.flight_starts_idx(i)-120:flightPaths.flight_starts_idx(i)+(10*120),:);
        day_index(i) = flightPaths.day(i);
        absTims(i) = flightPaths.flight_starts_idx(i);
    catch
        Flights2Use(:,:,i) =  flightPaths.AllFlights(flightPaths.flight_starts_idx(i-1)-120:flightPaths.flight_starts_idx(i-1)+(10*120),:);
        day_index(i) = flightPaths.day(i);
        absTims(i) = flightPaths.flight_starts_idx(i);
    end
end


% Get back to proper order:
Flights2Use = Flights2Use(:,:,out_markov.out_sort); % get the right flights in the right order
Flights2Use = Flights2Use(:,:,idx2use);
day_index = day_index(out_markov.out_sort);
day_index = day_index(idx2use);

absTims = absTims(out_markov.out_sort);
absTims = absTims(idx2use);



figure();
hold on;
c = jet(max(day_index));
for iii = 1:max(day_index)
    within_day =  find(day_index == iii);
    for ii = 1:length(within_day);
        ii = within_day(ii);
        plot(squeeze(Flights2Use(:,1,ii)),squeeze(Flights2Use(:,2,ii)),'color',c(iii,:));
    end
    
    clear within_day
end
colormap(jet)
colorbar();

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
    try
    tXYZ(i,1) = corr(Mflight(1:endOffset,1),squeeze(Flights2Use(1:endOffset,1,i)));
    tXYZ(i,2) = corr(Mflight(1:endOffset,2),squeeze(Flights2Use(1:endOffset,2,i)));
    tXYZ(i,3) = corr(Mflight(1:endOffset,3),squeeze(Flights2Use(1:endOffset,3,i)));
    catch
    tXYZ(i,1) = corr(Mflight(1:endOffset,1),squeeze(Flights2Use(1:endOffset,1,i-1)));
    tXYZ(i,2) = corr(Mflight(1:endOffset,2),squeeze(Flights2Use(1:endOffset,2,i-1)));
    tXYZ(i,3) = corr(Mflight(1:endOffset,3),squeeze(Flights2Use(1:endOffset,3,i-1)));
    end
end

% figure(); plot((tXYZ')');
% build smoothed line:
meav_vect = mean(tXYZ');
transitions2use = find(diff(day_index)>0);
transitions2use = [1 transitions2use];
for i = 1:length(transitions2use)-1;
    abs_new(transitions2use(i):transitions2use(i+1)) = absTims(transitions2use(i):transitions2use(i+1));
    gnew(transitions2use(i):transitions2use(i+1)) = smooth(meav_vect(transitions2use(i):transitions2use(i+1)),20);
    gnew(transitions2use(i)) = NaN;
end
    
  
figure();
hold on;
plot(mean(tXYZ')','.');
ylabel('Corr to median flight');
xlabel('flights (in order)')
% plot(smooth(mean(tXYZ'),15),'-','LineWidth',1);
% plot(gnew,'-','LineWidth',1);
% plot(1-mat2gray(diff(day_index)),'r')
    
figure();
hold on;
plot(mean(tXYZ')','.');
plot(gnew(isfinite(abs_new)),'-','LineWidth',1);
% plot(absTims,mean(tXYZ')','.');
% plot(abs_new(isfinite(abs_new)),gnew(isfinite(abs_new)),'-','LineWidth',1);
plot(1+(mat2gray(diff(day_index))/100),'--k')
ylabel('Corr to median flight');
xlabel('flights (in order)')

GG = mean(tXYZ')';
out.data = GG;
out.day_index = day_index;

% first and last quarter:
Q1 = mean(GG(1:round(size(GG,1)/4)));
Q4 = mean(GG((round(size(GG,1)/4))*2:round(size(GG,1)/4)*3-4));

