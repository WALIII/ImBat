for i = 1:length(ROIs_gal(1,:))
    figDir = dir(['*placeCell_' num2str(i) '.fig'])
    plotFiringTrajectoryAcrossDays(i) = figure('units','normalized','outerposition',[0 0 1 0.8]);
    sgtitle(['Gal ROI Selectivity 200311-200320: ROI ' num2str(i)]);
    for ii = 1:length(figDir)
        figure(plotFiringTrajectoryAcrossDays(i))
        h(ii)=subplot(ceil(length(ROIs_gal(:,ii))/3),3,ii);
        title(['Day ' num2str(ii)]);
        hold on;
        % Load saved figures
        figL(ii)=hgload(figDir(ii).name);
        % Prepare subplots
        % Paste figures on the subplots
        copyobj(allchild(get(figL(ii),'CurrentAxes')),h(ii));
        
    end
    saveas(plotFiringTrajectoryAcrossDays(i),[pwd filesep 'Gal_200311to20_placeCell_' num2str(i) '_acrossDays.tif']);
    savefig(plotFiringTrajectoryAcrossDays(i),[pwd filesep 'Gal_200311to20_placeCell_' num2str(i) '_acrossDays.fig']);
            
    close all;
end