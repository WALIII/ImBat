saveFlag = 1; %do you want to save the figures and output structure?
if saveFlag == 1
saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
% Check if folder exists
if exist([saveDir1 datestr(now,'yymmdd') filesep 'percSelectiveROI'])>0;
    disp('Youve been working today..');
else
    mkdir([saveDir1 datestr(now,'yymmdd') filesep 'percSelectiveROI'])
end
saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'percSelectiveROI' '\'];
end

plot_numSelectiveROI = figure(); 
subplot(1,2,1)
plot(squeeze(perc_flight_roi(1,1,:)))
hold on
plot(squeeze(perc_flight_roi(2,1,:)))
plot(squeeze(perc_flight_roi(3,1,:)))
legend('Pre','Dur','Post');
ylabel('% flights selective');
xlabel('Day');
set(gca,'xTick',[1:length(perc_flight_roi(1,1,:))]);
title('% of flights with at least one selective ROI');
hold off

subplot(1,2,2)
plot(squeeze(perc_flight_roi(1,2,:)))
hold on
plot(squeeze(perc_flight_roi(2,2,:)))
plot(squeeze(perc_flight_roi(3,2,:)))
legend('Pre','Dur','Post');
ylabel('% ROIs active');
xlabel('Day');
set(gca,'xTick',[1:length(perc_flight_roi(1,1,:))]);
title('% of unique selective ROIs');
hold off

if saveFlag == 1
    %save fig and tif of max projection
    %set(findall(maxFig,'-property','FontSize'),'FontSize',20);
    savefig(plot_numSelectiveROI,[saveDir 'scatterCorrDist_Gen_19to24_' datestr(now,'yymmdd-hhMMss') '.fig']);
    saveas(plot_numSelectiveROI, [saveDir 'scatterCorrDist_Gen_19to24_' datestr(now,'yymmdd-hhMMss') '.tif']);
    save([saveDir 'percSelectiveFlights_Gen_19to24_' datestr(now,'yymmdd-hhMMss') '.mat'],'perc_flight_roi');
end
