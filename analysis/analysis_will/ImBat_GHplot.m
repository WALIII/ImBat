function ImBat_GHplot(flightPaths);
% ImBat_GHplot.m

% GH Plot flights
% WAL3
% 2021.05


flights2plot = max(flightPaths.id);
cmap = jet(flights2plot);

figure();
hold on;
for i = 1:flights2plot;
    A = find(flightPaths.id == i);
    in = flightPaths.flight_ends_idx(A)-120; % 1s from flight end
    %in = flightPaths.flight_starts_idx(A); % plot flight starts

    plot(flightPaths.AllFlightsMasterTime([in;in]),[ones(size(in))+i;zeros(size(in))+i],'LineWidth',3,'color',cmap(i,:))
end

% plot reward:
for i = 1: size(flightPaths.RewardIdx,1);
    in = flightPaths.RewardIdx(i);
    in = flightPaths.AllFlightsMasterTime(in);
    plot([in;in],[ones(size(in));zeros(size(in))+flights2plot+1],'-r','LineWidth',0.03)
end
