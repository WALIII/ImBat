function [out] = ISI_Ephys(S_units,cell2use)



spikes_ts = {S_units.Timestamps};
xlim2use = 10000000;
% concat data together
idx = find([S_units.TT] == 4);

dat_combined = [];
for i = cell2use;
    clear spikes_ts1
spikes_ts1 = spikes_ts{idx(i)}/1000; % get spikes
spikes_ts1 = diff(spikes_ts1); % diff spikes
spikes_ts1(spikes_ts1>xlim2use) = []; % remove ISI above a threshold

dat_combined = cat(2,dat_combined,spikes_ts1); % combine data 
    edges_us = 10.^(0:0.05:6);
% 
% histogram((spikes_ts1),edges_us);
% set(gca,'XScale','log');
%     y1=get(gca,'ylim'); hold on; %plot([1e3 1e3],y1);  hold off;
%     xlabel('ms');   ylabel('ISI count');
%   %  title(['Average Firing Frequency: ' num2str(TT_unit(i).Frng_freq,3) ' Hz']);
% xlim([ 1000 xlim2use]);

%  set(gca, 'yscale','log')
%   set(gca, 'xscale','log')
end


% plot all
%figure(); 
histogram((dat_combined),edges_us,'edgecolor','none','FaceColor','black');
set(gca,'XScale','log');
    y1=get(gca,'ylim'); hold on; %plot([1e3 1e3],y1);  hold off;
    xlabel('ms');   ylabel('ISI count');
  %  title(['Average Firing Frequency: ' num2str(TT_unit(i).Frng_freq,3) ' Hz']);
% xlim([ 1000 xlim2use]);

out.dat = dat_combined;