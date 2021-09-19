function [TC_norm_pruned] = ImBat_MCS_helper_weight_ones(cs34_range,li)

    % Function to tabulate all 1s and calculate the non-false-absorbing state
    % thing
    flightpaths_in_TC_row = {}; 
    TC = zeros(li,li,2);
    for i=1:li
        if i==1
            continue
        else
            for m=1:length(cs34_range)
                if cs34_range(m) == i
                    one_ctr = 0;
                    for j=1:(length(cs34_range)-m)
                        if cs34_range(m+j) == 1
                            one_ctr = one_ctr + 1;
                        else
                            tc = cs34_range(m+j);
                            TC(i,tc,1) = TC(i,tc,1) + 1;
                            TC(i,tc,2) = TC(i,tc,2) + one_ctr;
                            break;
                        end
                    end
                end
            end
        end
    end

    TC_norm = zeros(li,li);
    for i=1:size(TC,1)
        row1 = TC(i,:,1);
        row2 = TC(i,:,2);
        row = row1./row2;
        for j=1:length(row)
            if isinf(row(j))
                row(j) = row1(j);
            end
        end
        row(isnan(row))=0;
        norm_row = row./(sum(row));
        TC_norm(i,:) = norm_row;
    end

    % Remove rows and columns w nan
    [TC_rows, TC_columns] = find(isnan(TC_norm));
    [TC_rows_nonnan, TC_columns_nonnan] = find(~isnan(TC_norm));
    flightpaths_in_TC_row{i} = unique(TC_rows_nonnan);
    TC_norm_pruned = TC_norm;
    TC_norm_pruned(unique(TC_rows),:) = [];
    TC_norm_pruned(:,unique(TC_rows)) = [];

    % Plot heatmap
    figure(); heatmap(TC_norm_pruned);
end
