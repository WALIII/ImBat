% Waterfall plots
% Load in cs45 

% Subset the c_s_34 data if needed
c_s_34_seg = c_s_34;
    
% Extract every 3-gram possible from the data using sliding window
TGM_all = [];
window_size=2;
for i=1:length(c_s_34_seg)-window_size
    three_gram = c_s_34_seg(i:i+window_size);
    TGM_all(i,:) = three_gram;
end

% Calculate co based on which flight types make up >5% of the data
five_percent = round(size(c_s_34,1)*(0.05));
co = 16;
colormap = distinguishable_colors(co);
length(c_s_34(c_s_34 > co));

% Remove all 3-grams that bridge sessions
TGM_blacklist = [];
[ss,rr] = sort(flightPaths34.flight_starts_idx);
day_idx = flightPaths34.day(rr);
dur_idx = flightPaths34.dur(rr);
for i=1:length(TGM_all)-2
    if day_idx(i) ~= day_idx(i+1)
        TGM_blacklist = [TGM_blacklist,[i-1,i,i+1]];
    end
end
TGM_pruned_1 = TGM_all;
TGM_pruned_1(TGM_blacklist,:) = [];

% Remove all 3-grams that contain flights greater than the cutoff 
TGM_pruned_2 = [];
for i=1:length(TGM_pruned_1)
    if sum(TGM_pruned_1(i,:) > co) >= 1
        continue
    else
        TGM_pruned_2 = [TGM_pruned_2;TGM_pruned_1(i,:)];
    end
end

% Remove all 3-grams containing a 1
TGM_pruned_3 = [];
for i=1:length(TGM_pruned_2)
    if sum(TGM_pruned_2(i,:)==1) >= 1 
        continue
    else
        TGM_pruned_3 = [TGM_pruned_3;TGM_pruned_2(i,:)];
    end
end

% Calculate mean duration of every flight type
% clustered_durs = {};
% for i=1:co
%     for m=1:length(c_s_34)
%         if c_s_34(m) == co
%             clustered_durs{i,:}= [clustered_durs(i),dur_idx(m)];
%         end
%     end
% end
            
% Make waterfall plots
for i=1:co
    TGM_bank = [];
    for j=1:length(TGM_pruned_2)
        if TGM_pruned_2(j,2) == i
            TGM_bank = [TGM_bank;TGM_pruned_2(j,:)];
        end
    end
    % Sort TGM bank by t-1 and t+1
    TGM_bank_sorted_1 = sortrows(TGM_bank,1);
    unique_tminus1 = unique(TGM_bank_sorted_1(:,1));
    TGM_sort3_temp = [];
    TGM_bank_sorted_2 = [];
    for j=1:length(unique_tminus1)
        sort_range = find(TGM_bank_sorted_1(:,1) == unique_tminus1(j));
        sort_tgm3 = sort(TGM_bank(sort_range(1):sort_range(end),3));
        TGM_sort3_temp = [TGM_sort3_temp;sort_tgm3];
    end
    TGM_bank_sorted_2 = [TGM_bank_sorted_1(:,1:2),TGM_sort3_temp];
    
    % Make the plot
    yaxis = [1:length(TGM_bank_sorted_2)];
    figure(); hold on; title(strcat("Waterfall plot cluster"," ",num2str(i)));
    xlabel("Trial (flight)"); ylabel("Time, aligned to takeoff of cluster")
    for j=1:length(TGM_bank_sorted_2)
        line([0 5],[yaxis(j) yaxis(j)], 'Color',colormap(TGM_bank_sorted_2(j,1),:),'LineWidth',4);
        line([5 10],[yaxis(j) yaxis(j)], 'Color',colormap(TGM_bank_sorted_2(j,2),:),'LineWidth',4);
        line([10 15],[yaxis(j) yaxis(j)], 'Color',colormap(TGM_bank_sorted_2(j,3),:),'LineWidth',4);
    end
    for j=1:length(TGM_bank_sorted_2)-1
        if TGM_bank_sorted_2(j,1) ~= TGM_bank_sorted_2(j+1,1)
            line([0 15],[yaxis(j) yaxis(j)], 'Color',[0.3 0.3 0.3],'LineWidth',2,'LineStyle',':');
        end
    end
    
    %pie plots
    pie_clusters = unique(TGM_bank_sorted_2(:,1));
    pie_bank = cell(length(pie_clusters),1);
    for j=1:length(unique(TGM_bank_sorted_2(:,1)))
        for m=1:length(TGM_bank_sorted_2)
            if TGM_bank_sorted_2(m,1) == pie_clusters(j)
                pie_bank{j} = [pie_bank{j};TGM_bank_sorted_2(m,:)];
            end
        end
    end
    
    for j=1:length(pie_bank)
        figure(); 
        pieData = pie_bank{j}(:,3);
        ax = gca(); 
        h = pie(ax, pieData); 
        % Define colors, one for each of the 3 wedges
        newColors = [colormap(pieData,:)];
        ax.Colormap = newColors; 
        delete(findobj(ax,'Type','text'))%option 2
        title(ax,strcat("Flights After Flight Type"," ",num2str(j),",",num2str(i)));
    end
end

figure(); hold on;
for i=1:co
    scatter(1,i,100,colormap(i,:),'filled');
end

