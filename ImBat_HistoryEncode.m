function [p_val p_val_pre] = ImBat_HistoryEncode(out);

% check stats for history encoding 


% run: [out] = ImBat_PlotMarkov(out_markov,FlightAlignedROI{1},48);
num_reps = 1000; % number of replications;
plotting_flag =0;


if plotting_flag ==1
figure(); 
tidx = (1:length(squeeze(out.Flightmat_post{5}(:,1,:))))/120;
for i = 1:3;
subplot(1,3,i)
hold on;
plot(tidx,squeeze(out.Flightmat_post{5}(:,i,:)),'b');
plot(tidx,squeeze(out.Flightmat_post{3}(:,i,:)),'r');
end
end

% get all flights w/ over 5 examples:
idx = cellfun('size',out.Flightmat_post,3);
idx2use = find(idx>4);
C = nchoosek(idx2use,2);

for i = 1:size(C,1);
%t = max(abs(mean(out.ROImat_post{C(i,1)}')-mean(out.ROImat_post{C(i,2)}')));
t = 1-corr(mean(out.ROImat_post{C(i,1)}')',mean(out.ROImat_post{C(i,2)}')');
clear combine_mat
combine_dat = cat(1,out.ROImat_post{C(i,1)}',out.ROImat_post{C(i,2)}');
for ii = 1:num_reps;
shuff_dat =  combine_dat(randperm(size(combine_dat,1)),:);
%combine_mat(ii) =  max(abs(mean(shuff_dat(1:size(out.ROImat_post{C(i,1)},2),:)-mean(shuff_dat(size(out.ROImat_post{C(i,1)},2):end,:)))));
combine_mat(ii) =  1-corr(mean(shuff_dat(1:size(out.ROImat_post{C(i,1)},2),:))', mean(shuff_dat(size(out.ROImat_post{C(i,1)},2):end,:))');
end

p = find(combine_mat>t);
p_val(i) = (size(p,2)/num_reps)*size(C,1); % for bonforoni correction
clear combine_mat combine_dat
end

% now look for pre-flight sig
% get all flights w/ over 5 examples:
idx2 = cellfun('size',out.Flightmat_pre,3);
idx2use2 = find(idx2>4);
C2 = nchoosek(idx2use2,2);

for i = 1:size(C2,1);
%t2 = max(abs(mean(out.ROImat_pre{C2(i,1)}')-mean(out.ROImat_pre{C2(i,2)}')));
t2 = 1-corr(mean(out.ROImat_pre{C2(i,1)}')',mean(out.ROImat_pre{C2(i,2)}')');

clear combine_mat
combine_dat2 = cat(1,out.ROImat_pre{C2(i,1)}',out.ROImat_pre{C2(i,2)}');
for ii = 1:num_reps;
shuff_dat2 =  combine_dat2(randperm(size(combine_dat2,1)),:);
combine_mat2(ii) =  max(abs(mean(shuff_dat2(1:size(out.ROImat_pre{C2(i,1)},2),:)-mean(shuff_dat2(size(out.ROImat_pre{C2(i,1)},2):end,:)))));
%combine_mat2(ii) =  max(abs(mean(shuff_dat2(1:size(out.ROImat_pre{C2(i,1)},2),:)-mean(shuff_dat2(size(out.ROImat_pre{C2(i,1)},2):end,:)))));
end

p2 = find(combine_mat2>t2);
p_val_pre(i) = (size(p2,2)/num_reps)*size(C2,1); % for bonforoni correction
clear combine_mat combine_dat
end