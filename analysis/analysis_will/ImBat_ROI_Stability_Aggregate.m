function ImBat_ROI_Stability_Aggregate

% Temp function to aggregate and quantify ROI stability over days

% WAL3
% 02/23/2021

% nb: Run in data in dir with 'processed' folders...


%% Aggregate Data

DIR = pwd;
files = dir(pwd);
files(ismember( {files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.

for i = 1 : length(subFolders)
    load([subFolders(i).name, '/CellReg_files/ROI_Data/Saved_Data/Aligned_Data.mat'],'CombinedROI','flightPaths')
    %load([subFolders(i).name, '/CellReg_files/ROI_Data/ROI_Data.mat'])
    [out{i}] = ImBat_analysis_11122020(CombinedROI,flightPaths);
    clear CombinedROI flightPaths
   disp('Moving to next session...'); 
end

% Initialize vars:
out_final.x_2save = [];
out_final.y_2save = [];
out_final.err_2save = [];
out_final.FL = [];

out_final2.x_2save = [];
out_final2.y_2save = [];
out_final2.err_2save = [];
out_final2.FL = [];
 
%% Concat data
for ii = 1:size(out,2);
   out_final.x_2save = cat(1,out_final.x_2save,out{ii}.forward.x_2save);
   out_final.y_2save = cat(1,out_final.y_2save,out{ii}.forward.x_2save);
   out_final.err_2save = cat(1,out_final.err_2save,out{ii}.forward.err_2save);
   out_final.FL = cat(1,out_final.FL,out{ii}.forward.FL);
% Reverse Sort
   out_final2.x_2save = cat(1,out_final2.x_2save,out{ii}.reverse.x_2save);
   out_final2.y_2save = cat(1,out_final2.y_2save,out{ii}.reverse.x_2save);
   out_final2.err_2save = cat(1,out_final2.err_2save,out{ii}.reverse.err_2save);
   out_final2.FL = cat(1,out_final2.FL,out{ii}.reverse.FL);

end
        out_final.ROI_ON = out{1}.forward.ROI_ON;
        out_final2.ROI_ON = out{1}.reverse.ROI_ON;

% Sort data:
[~, sort_idx] = sort(out_final.x_2save(:,1),'ascend');
out_final.x_2save = out_final.x_2save(sort_idx,:);
out_final.y_2save = out_final.y_2save (sort_idx,:);
out_final.err_2save = out_final.err_2save(sort_idx,:);
out_final.FL = out_final.FL(sort_idx,:);

[~, sort_idx2] = sort(out_final2.x_2save(:,1),'descend');
out_final2.x_2save = out_final2.x_2save(sort_idx2,:);
out_final2.y_2save = out_final2.y_2save (sort_idx2,:);
out_final2.err_2save = out_final2.err_2save(sort_idx2,:);
out_final2.FL = out_final2.FL(sort_idx2,:);
        


% make plot:
ImBat_PlotRedBlue(out_final)
ImBat_PlotRedBlue(out_final2)
xlim([-200 300])


%% remove high and low vals( for stats)

out_final3 = out_final;
out_final4 = out_final2;

dat2remove = out_final3.x_2save(:,1)-out_final3.ROI_ON<-30 | out_final3.x_2save(:,1)-out_final3.ROI_ON-(out_final3.FL+60)>0;

out_final3.y_2save(dat2remove,:) = [];
out_final3.err_2save(dat2remove,:) = [];
out_final3.FL(dat2remove,:) = [];
out_final3.x_2save(dat2remove,:) = [];

out_final4.y_2save(dat2remove,:) = [];
out_final4.err_2save(dat2remove,:) = [];
out_final4.FL(dat2remove,:) = [];
out_final4.x_2save(dat2remove,:) = [];

ImBat_PlotRedBlue(out_final3);
title('removing pre and post ROIs');

D2  = ImBat_PlotRedBlue(out_final4);

% % now reverse the order:
% % Now, align to reward, and reverse sort:
% out_final4 = out_final3;
% offset2use = 120;
% out_final4.ROI_ON = out_final4.ROI_ON+offset2use;;
% out_final4.x_2save = out_final4.x_2save-out_final4.FL+offset2use-20;
% [~, sort_idx] = sort(out_final4.x_2save(:,1),'descend');
% out_final4.x_2save = out_final4.x_2save(sort_idx,:);
% out_final4.y_2save = out_final4.y_2save (sort_idx,:);
% out_final4.err_2save = out_final4.err_2save(sort_idx,:);
% out_final4.FL = out_final4.FL(sort_idx,:);
% D2 = ImBat_PlotRedBlue(out_final4);
% title('removing pre and post ROIs');


%% Shuffle and more stats:
out_final5 = out_final;
sort_idx = randperm(size(out_final5.x_2save,1));
out_final5.x_2save(:,2) = out_final5.x_2save(sort_idx,2);
D1 = ImBat_PlotRedBlue(out_final5);

% figure(); 

% % Jitter Field
% out_final6 = out_final;
% out_final6.x_2save(:,2) = out_final6.x_2save(sort_idx,2)-(randi(10*30,size(out_final6.x_2save,1),1)-5*30);
% D3 = ImBat_PlotRedBlue(out_final6);



 E1 = D1.X1;
 E2 = D2.X1;
 E1(E1>8) = [];
 E2(E2>8) = [];
 
 figure();
 hold on;
histogram((E1),30,'Normalization','probability')
% histogram(mat2gray(D3.X1),50,'Normalization','probability')
histogram((E2),30,'Normalization','probability')
title('Timing change in field location');
ylabel('Incidence');
xlabel('Time (s)')
legend('Shuffle', 'Longest Tracked'); 



