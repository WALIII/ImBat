function [out] = ISI_Calcium(CombinedROI,cell2use,day2use);

% Calculate and plot the inter-spike-interval (isi) for the calcium data

% user inputs: 
%cell2use = 1;%size(CombinedROI.S,1)';
% day2use = 1;

Spk = 0;
B33 = 0;
for i = cell2use
    % binarie spikes
    ind2use = find(CombinedROI.day_vector ==day2use);
    Spikes = CombinedROI.S(i,ind2use);
    Spikes = zscore(Spikes);
    thresh = 3;

    Spikes(Spikes>thresh) = thresh;
   ind2rm = find(Spikes<thresh);
    Spk_temp = find(Spikes ==thresh);
    Spk = [Spk max(Spk)+Spk_temp];
    clear Spikes
    % get ISI
    Spikes = CombinedROI.S(i,:);
    Spikes(ind2rm) = 0;
    scale_Spikes= Spikes*1/min(Spikes(Spk_temp));
    below_33ms = sum(scale_Spikes)/2;
    B33 = B33+below_33ms;
    clear Spikes ind2rm below_33ms scale_Spikes
end

Spk = smooth(Spk,3)';
x = diff(Spk./(1/30));
 x = [x ones(1,round(B33/10))*33];
x_rand = randi(30,1,size(x,2))-10;
x = x+x_rand;
x(x<35) = [];

%mean(1./x)*100
% % 
% % figure(); 
% x = diff(Spk);
% x = [x ones(1,round(B33))*.15];
% x(x>60) = [];
% histogram(x,'BinWidth',1,'Normalization','probability');
% set(gca,'YScale','log')
% %set(gca, 'xscale','log')
% 
% xt = xticklabels;
% for i = 1:7;
% xt{i} = num2str((i-1)*10 *33);
% end
% xticklabels(xt)




% figure();
% histogram(x,'BinWidth',1,'Normalization','probability');
% set(gca,'XScale','log')

edges_us = 10.^(0:0.05:6);

histogram((x),edges_us,'edgecolor','none','FaceColor','black');
set(gca,'XScale','log')
   % xlabel('ms');  % ylabel('ISI count');

 out.dat = x;
 out.fr = mean(1./x)*100;
