function out = ImBat_FlightStats(FlightAlignedROI);

% Get fight stas over time for the top 3-4 flights..
plot_data = 1;

% Compute PCA:
t2use = ['Sorting based on PCA'];
    VV_full = squeeze(FlightAlignedROI.ClustFlight(50:550,:,:));

    %% PCA!!
    VV_full2 = permute(VV_full,[1,3,2]); % permute order
    F2PCA = reshape(VV_full2,size(VV_full,1)*size(VV_full2,2),size(VV_full2,3)); % reshape matrix
    [coeff,score,latent,tsquared,explained,mu] = pca(F2PCA); % do PCA
    VV = reshape(score(:,1),size(VV_full,1),size(VV_full,3)); % reshape matrix
    
    explained(1) % report % explained...
    
    
  dayindex = 1;  
  FirstDay = 0;
  move = 0;
for i = 1:max(FlightAlignedROI.CutCells_date)
    cF = find(FlightAlignedROI.CutCells_date==i);
    if FirstDay==0 &&  size(cF,2)>0;
        FirstDay = i;
        move =1;
    elseif move ==0;
        continue
    end
    tdF = squeeze(VV(:,cF));
    if i ==FirstDay;
    meanFlight = mean(tdF,2);
    end
    for ii = 1: size(tdF,2);
        R{dayindex}(ii) = corr2(meanFlight,tdF(:,ii));
    end
    
    try
    toPlot(1,dayindex) = mean(R{dayindex});
    toPlot(2,dayindex) = std(R{dayindex})/2;
    toPlot(3,dayindex) = size(tdF,2);
    clear tdF cF
    [pval_combined_data(i),~] = ranksum((R{1}(:)), (R{dayindex}(:)),'tail','right');
    catch
         toPlot(1,dayindex) = nan;
    toPlot(2,dayindex) = nan;
    pval_combined_data(dayindex) = 1;
    end
    dayindex = dayindex+1;
end


out.toPlot = toPlot;
out.pval_combined_data = pval_combined_data;

if plot_data ==1;
    figure();
    col = 'r';
hold on;
 plot(toPlot(1,:),col);

    y1 = toPlot(1,:);
    x1 = 1:length(y1);
    err1 = toPlot(2,:);
errorbar(x1,y1,err1,'vertical','o','color',col,'CapSize',1);

% figure(); plot(pval_combined_data);

[a b] = find(pval_combined_data<(0.05./i));
scatter(b,ones(length(b),1),'k*')
title('flight correlation vs time');
ylabel(' Correlation to day 1');
xlabel(' days');
else
    disp('Plotting turned off');
end

