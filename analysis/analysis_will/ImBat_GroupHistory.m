function ImBat_GroupHistory
% ImBat_GroupHistory

% Run in 'Processed' folder.


% index into the right folder
Mf = [];
Mn = [];
Mf2 = [];
Mn2 = [];
DIR = pwd;
Flights_Mean_perDay_Temp = [];

files = dir(pwd);
files(ismember( {files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.


for i = 1 : length(subFolders)
        load([subFolders(i).name, '/CellReg_files/ROI_Data/Saved_Data/Aligned_Data.mat'])

        % calculate transition probs:
        [out_markov] = ImBat_New_Markov(flightPaths);
        close all

 [Mf_temp, Mn_temp, Mf2_temp, Mn2_temp,out_dat{i}] = history_scrap2(out_markov,FlightAlignedROI);
 Mf = cat(2,Mf, Mf_temp); %sum
 Mn = cat(2,Mn, Mn_temp);
 Mf2 = cat(2,Mf2, Mf2_temp); %size
 Mn2 = cat(2,Mn2, Mn2_temp);

 
end

% plotting; 
figure();
hold on;
for i = 1:size(Mf,2);
x = [ 1 2];
y = [Mf(i)./Mf2(i) Mn(i)./Mn2(i)];;
plot(x, y, '-k')
hold on
scatter(x(1), y(1), 50, 'b', 'filled')
scatter(x(2), y(2), 50, 'r', 'filled')
end


end
% subfuncitons 

function [Mf, Mn, Mf2, Mn2, out_dat] = history_scrap2(out_markov,FlightAlignedROI)
alphaVal = 0.05;

for iii = 1:2;
clear sig_cell_pre  sig_cell_post idx
counter = 1;
for i = 1:size(FlightAlignedROI{1}.C,1);
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

try
    Mf(iii) = sum(sig_cell_pre);
    Mn(iii) = sum(sig_cell_post);
    
    Mf2(iii) = size(sig_cell_pre,2);
    Mn2(iii) = size(sig_cell_post,2);
    
    % export data
    out_dat{iii}.sig_cell_pre = sig_cell_pre;
    out_dat{iii}.sig_cell_post = sig_cell_post;
    out_dat{iii}.idx = idx;
catch
    Mf(iii) = NaN;
    Mn(iii) = NaN;
    
    Mf2(iii) = NaN;
    Mn2(iii) = NaN;
    
    % export data
    out_dat{iii}.sig_cell_pre = NaN;
    out_dat{iii}.sig_cell_post = NaN;
    out_dat{iii}.idx = NaN;
end
end
end



% Plotting


 