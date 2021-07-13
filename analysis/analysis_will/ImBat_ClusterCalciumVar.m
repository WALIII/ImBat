function [pval_combined_data_DR,c, c_sort]= ImBat_ClusterCalciumVar(FlightAlignedROI,roi2plot);

CutCells = FlightAlignedROI.C_raw;
CutFlights = FlightAlignedROI.ClustFlight;
cell2plot = zscore(squeeze(CutCells(roi2plot,:,:)));
toRemove = sum(cell2plot);
ind2Remove = find(toRemove==0);
cell2plot(:,ind2Remove) = [];
CutFlights(:,:,ind2Remove) = [];
ROI_ON  = 150;
bound2plot = 1:350;
clust2use = 3;
pval_combined_data_DR = 1;

ROI2Compare2 = 1;
try
l = linkage(cell2plot(150:250,:)','ward');
catch
    disp('NO DATA');
    return
end
c=cluster(l,'maxclust',clust2use);
[aa,bb]=sort(c);
bound2use = diff(aa);
line2plot = find(bound2use==1);
figure(); hold on;
imagesc(cell2plot(bound2plot,bb)');
colormap(parula.^2);
  ax = gca;
  c_sort = bb; % the sort 
    % Set where ticks will be
    ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
    % Set TickLabels;
    ylabel('filghts')
    xlabel('time from takeoff');
    ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};
    
for i = 1:size(line2plot,1)
line([1 max(bound2plot)],[line2plot(i) line2plot(i)],'color','r','LineWidth',3)
end
axis tight
figure(); 
for i = 1:clust2use
subplot(2,3,i);
ind2use = find(c==i);
imagesc(cell2plot(bound2plot,ind2use)');
  ax = gca;
    % Set where ticks will be
    ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
    % Set TickLabels;
    ylabel('filghts')
    xlabel('time from takeoff');
    ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};
end
colormap(parula.^2);



% now plot the flights, against shuffled data
figure();
% subplot(2,3,4:6);
hold on;
col = {[1, 0, 0, 0.2],[0, 1, 0, 0.2],[0, 0, 1, 0.2],'c','m','k'};
for i = 1:clust2use
ind2use = find(c==i);
fakeInd = randperm((size(CutFlights,3)));
%fakeInd = fakeInd(1:size(ind2use,1)*2);
plot3(squeeze(CutFlights(:,1,ind2use)),squeeze(CutFlights(:,2,ind2use)),squeeze(CutFlights(:,3,ind2use)),'color',col{i},'LineWidth',2);
A{i}(:,:,:) = squeeze(CutFlights(1:550,:,ind2use));
Afake{i}(:,:,:) = squeeze(CutFlights(1:550,:,fakeInd));


end

Afake2{i}(:,:,:) = squeeze(CutFlights(1:550,:,fakeInd));

% get Euclidian distance between each group and its mean,
for i =ROI2Compare2;
Mtrue = median(squeeze(A{i}(:,:,:)),3);
Mshuff = mean(squeeze(Afake{i}(:,:,:)),3);
   for ii = 1:size(A{i},3)
            % Euclidian distance
            FL(ii) = sum(sqrt(sum(squeeze(A{i}(:,:,ii))-Mtrue(:,:)) .^ 2));
            FL_shuff(ii) = sum(sqrt(sum(squeeze(Afake{i}(:,:,ii))-Mshuff(:,:)).^ 2));
   end
   for ii = 1:size(A{clust2use},3) % change this to 2, o3 for the histogram 
            % Euclidian distance
               FL2(ii) = sum(sqrt(sum(squeeze(A{clust2use}(:,:,ii))-Mtrue(:,:)) .^ 2));
   end
%      for ii = 1:size(A{1},3)
%             % Euclidian distance
%                FL3(ii) = sum(sqrt(sum(squeeze(A{1}(:,:,ii))-Mtrue(:,:)) .^ 2));
%      end
end

figure();
hold on;
histogram(FL,20,'Normalization','probability','BinWidth',10000,'FaceColor','b');
histogram(FL2,20,'Normalization','probability','BinWidth',10000,'FaceColor','g');
%histogram(FL3,20,'Normalization','probability','BinWidth',10000,'FaceColor','r');

[pval_combined_data_DR,~] = ranksum(FL,FL2)
% All-all corr



% figure(); 
% hold on;
% D  = sqrt(sum((squeeze(A(:,1,:))' - squeeze(A(:,3,:))') .^ 2));
% D2  = sqrt(sum((squeeze(A(:,1,:))' - squeeze(A(:,2,:))') .^ 2));
% D3  = sqrt(sum((squeeze(A(:,3,:))' - squeeze(A(:,2,:))') .^ 2));
% % Shuffled data 
% fD  = sqrt(sum((squeeze(Afake(:,1,:))' - squeeze(Afake(:,3,:))') .^ 2));
% fD2  = sqrt(sum((squeeze(Afake(:,1,:))' - squeeze(Afake(:,2,:))') .^ 2));
% fD3  = sqrt(sum((squeeze(Afake(:,3,:))' - squeeze(Afake(:,2,:))') .^ 2));
% 
% MD = (D+D2+D3)/3;
% fakeMD = (fD+fD2+fD3)/3;
% plot(MD./fakeMD);
% 
% figure(); 
% hold on;
% plot(MD); 
% plot(fakeMD);
% 
% % make a distribution of fakes:
% 
% for ii = 1: 1000;
% for i = 1:3
% ind2use = find(c==i);
% fakeInd = randi(size(CutFlights,3),size(ind2use,1),1);
% A(:,i,1) = mean(squeeze(CutFlights(:,1,ind2use)),2);
% A(:,i,2) = mean(squeeze(CutFlights(:,2,ind2use)),2);
% A(:,i,3) = mean(squeeze(CutFlights(:,3,ind2use)),2);
% 
% Afake(:,i,1) = mean(squeeze(CutFlights(:,1,fakeInd)),2);
% Afake(:,i,2) = mean(squeeze(CutFlights(:,2,fakeInd)),2);
% Afake(:,i,3) = mean(squeeze(CutFlights(:,3,fakeInd)),2);
% 
% end
% 
% D  = sqrt(sum((squeeze(A(:,1,:))' - squeeze(A(:,3,:))') .^ 2));
% D2  = sqrt(sum((squeeze(A(:,1,:))' - squeeze(A(:,2,:))') .^ 2));
% D3  = sqrt(sum((squeeze(A(:,3,:))' - squeeze(A(:,2,:))') .^ 2));
% % Shuffled data 
% fD  = sqrt(sum((squeeze(Afake(:,1,:))' - squeeze(Afake(:,3,:))') .^ 2));
% fD2  = sqrt(sum((squeeze(Afake(:,1,:))' - squeeze(Afake(:,2,:))') .^ 2));
% fD3  = sqrt(sum((squeeze(Afake(:,3,:))' - squeeze(Afake(:,2,:))') .^ 2));
% 
% MD = (D+D2+D3)/3;
% fakeMD = (fD+fD2+fD3)/3;
% fakeDist(ii,:)= (MD./fakeMD);
% end
% 
% col = jet(2);
% figure();
% hold on;
% adata = fakeDist;
% L = size(adata,2);
% se = std(adata)/4;%sqrt(length(adata));
% mn = mean(adata);
% mn = smooth(mn,1)';
% h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col(1,:)); alpha(0.5);
% plot(mn,'Color',col(1,:));
% 
% plot(zscore(A(:,1,1)),'Color',col(1,:));
% 
% figure(); plot(diff(mn));
% 
% 
% 
% 
% 
% 
