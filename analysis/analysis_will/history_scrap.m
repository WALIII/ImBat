function history_scrap(out_markov,FlightAlignedROI)
alphaVal = 0.05;

for iii = 1:2;
clear sig_cell_pre  sig_cell_post idx
counter = 1;
for i = 1:50;
    close all;
    clear out p_val_pre p_val;
    
    
%     sig_cell_pre(counter) = 0; % future/ planning
%     sig_cell_post(counter) = 0; %history
    
    [out] = ImBat_PlotMarkov(out_markov,FlightAlignedROI{iii},i);
    try
        [p_val, p_val_pre] = ImBat_HistoryEncode(out);
        
        % store data for later
        store_data.pval{i} = p_val;
        store_data.pval{i} = p_val;

        if size(find(p_val<alphaVal),2)>0;
            sig_cell_post(counter) = 1;
        else
            sig_cell_post(counter) = 0;
        end
        
        if size(find(p_val_pre<alphaVal),2)>0;
            sig_cell_pre(counter) = 1;
        else
            sig_cell_pre(counter) = 0;
        end
        idx(counter) = i;
        counter = counter+ 1;
        
    catch
        disp('Not enough data to compare');
    end
    
end

    Mf(iii) = sum(sig_cell_pre);
    Mn(iii) = sum(sig_cell_post);
    
    Mf2(iii) = size(sig_cell_pre,2);
    Mn2(iii) = size(sig_cell_post,2);
    
    % export data
    out_dat{iii}.sig_cell_pre = sig_cell_pre;
    out_dat{iii}.sig_cell_post = sig_cell_post;
    out_dat{iii}.idx = idx;
    
    end
% Plotting


 
figure();
hold on;
for i = 1:size(Mf2,2);
x = [ 1 2];
y = [Mf(i) Mn(i)];;
plot(x, y, '-k')
hold on
scatter(x(1), y(1), 50, 'b', 'filled')
scatter(x(2), y(2), 50, 'r', 'filled')
 
end