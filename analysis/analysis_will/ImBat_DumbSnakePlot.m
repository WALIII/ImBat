function ImBat_DumbSnakePlot();

% run in processed folder


DIR = pwd;

files = dir(pwd);
files(ismember( {files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.

disp('getting data...')
for i = 1 : length(subFolders)

load([subFolders(i).name, '/CellReg_files/ROI_Data/Saved_Data/Aligned_Data.mat'])

    
    % get all 'FlightAlignedROI' data
    FL_Struct{i} = FlightAlignedROI;
    
    clear FlightAlignedROI
end


% Plot the damn data

G = [];
for iii = 1:length(FL_Struct)
    clear FlightAlignedROI;
FlightAlignedROI = FL_Struct{iii};
for i = 1:3
    GG = FlightAlignedROI{i}.S;
    G = cat(1,G,zscore(squeeze(mean(GG,3)'))');
    clear GG
end

end

% remove zeros
Gm = sum(G');
ind2rmv = find(Gm ==0); %remove zeros...
    G(ind2rmv,:) = []; 
    
    


    G_smth = [];
    
% smooth data
smth_order = 20;
for ix = 1:size(G,1);
    G_smth(ix,:) = smooth(G(ix,:),smth_order);
end

G_smth = zscore(G_smth')';
[~, gm] = max(G_smth(:,50:350),[],2);
[gms idx] = sort(gm,'descend');
% plot for real

ROI_ON = 125;
figure();
hold on;
imagesc(G_smth(idx,:),[0 6]); 
colormap(hot);
plot([ROI_ON ROI_ON],[0 size(G_smth,1)],'--w')
axis tight
xlim([0 450])
ax = gca;

ax.XTick = [ ROI_ON-90 ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60 ROI_ON+90 ROI_ON+120,ROI_ON+150,ROI_ON+180,ROI_ON+210,ROI_ON+240];
% Set TickLabels;
ylabel('trajectory aligned ROI activity')
xlabel('Time from takeoff');
ax.XTickLabel = {'-3','-2','-1','0','1','2','3','4','5','6','7','8'};

colorbar();



