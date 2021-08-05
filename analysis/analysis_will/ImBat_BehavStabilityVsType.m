
function ImBat_BehavStabilityVsType(out_markov,flightPaths,flight2use,toPlot)
% look at the stability of behavior over time


% colormaps
cmap3 = colormap(lines(7));
cmap1 = colormap(cubehelix((max(out_markov.FlightIDVector)-7),1.5,3,4,1,[0.29,0.92]));

col = cat(1,cmap3,cmap1);
col = cat(1,[0.7 0.7 0.7],col);


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


pre_sort_idx = out_markov.FlightIDVector(idx2use+1);
%[a2 b2] = sort(out_markov.FlightIDVector(idx2use-1),'ascend');
% idx2use = idx2use(b2);

G1 = unique(pre_sort_idx);

if exist('toPlot')
    G1 = toPlot;
else
    toPlot = 1:size(G1,2);
end

plot_all_flights = 0;
figure(); hold on;
for i = 1:size(G1,2);
    id_temp = find(pre_sort_idx== G1(i));
    Mflight = mean(Flights2Use(:,:,id_temp),3); % mean flight
    %plot(Mflight(:,3),'color',col(G1(i),:),'LineWidth',2);
    plot3(Mflight(:,1),Mflight(:,2),Mflight(:,3),'color',col(G1(i),:),'LineWidth',2);
    if plot_all_flights ==1;
        for ii = 1:size(Flights2Use(:,:,id_temp),3)
            plot3(squeeze(Flights2Use(:,1,id_temp)),squeeze(Flights2Use(:,2,id_temp)),squeeze(Flights2Use(:,3,id_temp)),'color',col(G1(i),:),'LineWidth',.5);
        end
    end
    clear id_temp
    
end




figure();
for ii = 1:3;
    subplot(1,3,ii);
    hold on;
    for i = 1:size(G1,2);
        try
            id_temp = find(pre_sort_idx== G1(i));
            
            adata = squeeze(Flights2Use(:,ii,id_temp))'; % mean flight
            L = size(adata,2);
            se = nanstd(adata);%sqrt(length(adata));
            mn = nanmedian(adata);
            mn = smooth(mn,10)';
            se = smooth(se,10)';
            h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col(G1(i),:)); alpha(0.5);
            plot(mn,'Color',col(G1(i),:));
            clear id_temp adata
        catch
        end
    end
end




% figure(); hold on;
% for i = [2 3 5];
%     try
%    id_temp = find(pre_sort_idx== G1(i));
%
%     adata = squeeze(Flights2Use(:,3,id_temp))'; % mean flight
% L = size(adata,2);
% se = nanstd(adata)/2;%sqrt(length(adata));
% mn = nanmedian(adata);
% mn = smooth(mn,10)';
% se = smooth(se,10)';
% h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col(G1(i),:)); alpha(0.5);
% plot(mn,'Color',col(G1(i),:));
% clear id_temp adata
%     catch;
%     end
% end
