function [Lcorr, McorrM,Rcorr] = Scrap_remap(out);


% out = ImBat_PlotAlignedROIs(FlightAlignedROI{1},ROI_Data,flightPaths,'sort',dark_cluster{2},'smooth',1); % cluster 2 is # 2 in FLightAlighedROI

minNumFlights = 5; % min num of flights to compare
counter = 1;

for i = 1:size(out.ROI_sorted,2);
clear Light Dark
    % get transition point
A1 = out.ROI_sorted{i};

% get Light/Dark
Dark = A1(1:out.transition(2),:);
Light = A1(out.transition:end,:);

% Dark = A1(1:out.transition(2),1:350);
% Light = A1(out.transition:end,1:350);
% Remap = Rem1(out.transition:end,1:350);

% remove zeros

Dark(find(mean(Dark')==0),:) = [];
Light(find(mean(Light')==0),:) = [];


% get remapped data:
Remap = [];
 while size(Remap,1)<minNumFlights
     disp('Remap too small, resampling...');
Rem1 = out.ROI_sorted{randi(size(out.ROI_sorted,2),1)};
Remap = Rem1(out.transition:end,:);
Remap(find(mean(Remap')==0),:) = [];
 end



%if there are a min number in each group:
if size(Dark,1)>minNumFlights & size(Light,1)>minNumFlights*2 

 mDark =  mean(Dark);
 
 mLighto = mean(Light(1:2:end,:));
 mLighte = mean(Light(2:2:end,:));
 mRemape = mean(Remap(2:2:end,:));
 
 % calculate corr:
 temp = corrcoef(mLighto,mLighte);
 Lcorr(counter) = temp(2);
 temp = corrcoef(mLighto,mDark);
 Mcorro(counter) = temp(2);
  temp = corrcoef(mLighte,mDark);
 Mcorre(counter) = temp(2);
McorrM(counter) = (Mcorro(counter)+Mcorre(counter))/2;
% remap
 temp = corrcoef(mRemape,mLighte);
 Rcorr(counter) = temp(2);



 counter = counter+1;
end


 
end

try
figure(); 
hold on;
histogram(1-McorrM,'BinWidth',0.15,'Normalization','probability','FaceColor','k');
histogram(1-Lcorr,'BinWidth',0.15,'Normalization','probability','FaceColor','c');
histogram(1-Rcorr,'BinWidth',0.15,'Normalization','probability','FaceColor','r');

[pval_combined_data,~] = ranksum(McorrM,Lcorr)
catch
    disp( ' Not enough trials');
    McorrM = 1;
    Lcorr = 1;
    Rcorr =1;
end