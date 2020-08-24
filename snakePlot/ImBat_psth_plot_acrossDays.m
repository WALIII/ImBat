function ImBat_psth_plot_acrossDays(psth_stableROI)
saveFlag = 1;
plotDir = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
colDay = jet(length(psth_stableROI.velocity));
% Check if folder exists
if exist([plotDir datestr(now,'yymmdd')])>0;
    disp('Youve been working today..');
else
    mkdir([plotDir datestr(now,'yymmdd')])
end
saveDir = [plotDir datestr(now,'yymmdd') '\'];
%cd(saveDir);
plotPSTH_acrossDays = figure('units','normalized','outerposition',[0 0.1 0.4 0.8]);
for cell_ii = 1:size(psth_stableROI.traceMean{1}{1},1)
    for day_ii = 1:length(psth_stableROI.velocity)
     %try
            %plot flight paths of each cluster
            nFlights = size(psth_stableROI.velocity{day_ii}{2}{1},1); %number of flights in this cluster
            colFlight = jet(nFlights); %color each cluster differently)
            for flight_i = 1:nFlights
                plotPaths = subplot(10,1,[1 2]);
                plot3(squeeze(psth_stableROI.flight{day_ii}{2}{cell_ii}(flight_i,1,:)),squeeze(psth_stableROI.flight{day_ii}{2}{cell_ii}(flight_i,2,:)),squeeze(psth_stableROI.flight{day_ii}{2}{cell_ii}(flight_i,3,:)),...
                    '-','LineWidth',1,'Color',colDay(day_ii,:));
                hold on;
                view(0,90);
                xlim([-3 3])
                ylim([-3 3])
                title(['ROI # ' num2str(cell_ii) ': Gal 200311to20 - clust2']);
                    xlabel('m'); ylabel('m');
                
                %plot velocity of each cluster
                plotVelocity = subplot(10,1,[3 4]);
                plot(1:length(psth_stableROI.velocity{day_ii}{2}{cell_ii}(flight_i,:)),psth_stableROI.velocity{day_ii}{2}{cell_ii}(flight_i,:),'Color',colDay(day_ii,:));
                hold on;
                    ylabel('Velocity (m/s)');
                    yt = get(gca,'YTick');
                    xt = get(gca,'XTick');
                    set(gca,'xticklabel',[]);
                    %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
           
                ylim([0 4.5]);
                xlim([0 length(psth_stableROI.velocity{1}{2}{cell_ii}(1,:))]);
                
                end
            %plot mean and stdev of activity per cluster
            plotMeanTraces = subplot(10,1,[5 10]);
            boundedline(1:length(psth_stableROI.traceMean{day_ii}{2}(cell_ii,:)),psth_stableROI.traceMean{day_ii}{2}(cell_ii,:),psth_stableROI.traceStd{day_ii}{2}(cell_ii,:),...
                'alpha','cmap',colDay(day_ii,:),'nan','gap');
            hold on;
            ylabel('Mean stdev df/f');
            set(gca,'xticklabel',[]);
            xlim([0 length(psth_stableROI.traceMean{1}{2}(cell_ii,:))]);
            xlabel('Time (s)');
                xt = get(gca, 'XTick');
                set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
       % catch
           % sgtitle(['Flight Aligned Activity per cluster: Gal All Cell# ' num2str(cell_ii) ' NO ACTIVITY']);
       % end   
    
    end
    drawnow;
    if saveFlag == 1
    saveas(plotPSTH_acrossDays,[saveDir filesep 'Gal_200311to20_psth_acrossDays_ROI' num2str(cell_ii) datestr(now,'yymmdd_HHMM') '.tif']);
   savefig(plotPSTH_acrossDays,[saveDir filesep 'Gal_200311to20_psth_acrossDays_ROI' num2str(cell_ii) datestr(now,'yymmdd_HHMM') '.fig']);
    end
   clf;
end