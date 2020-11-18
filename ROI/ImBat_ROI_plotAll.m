function ImBat_ROI_plotAll
close all;
batId = 'Z';
saveFlag = 1;
%make saving directory
if saveFlag == 1
    saveTag = '';
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    %saveDir1 = '/Users/periscope/Desktop/analysis/flight/plots/';
    %saveDir1 = 'C:\Users\tobias\Desktop\analysis\plots\';
    if ~exist([saveDir1 datestr(now,'yymmdd') filesep 'maxProj_allDays'])
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'maxProj_allDays']);
    else
        disp('You have been working today...');
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'maxProj_allDays' filesep];
end


dirTop = dir([batId '*']);
for bn_i = 1:length(dirTop)
    batInit(bn_i,1:2) = dirTop(bn_i).name(1:2);
end
uniqueInit = unique(batInit);

for bat_i = 1:length(uniqueInit)-1
    close all;
    dirBat = dir([batId uniqueInit(bat_i) '*']);%uniqueInit(bat_i+1) '*']);
    plotAllROI =  figure('units','normalized','outerposition',[0 0 1 1]);
    sgtitle([dirBat(1).name  ': Max Projections Full Session']);
    %ha = tight_subplot(ceil(length(dirBat)/5),5,[.02 .01],[.01 .08],[.01 .01]);
    for day_i = 1:length(dirBat);
        cd(dirBat(day_i).name);
        dirFly = dir('*fly*extraction');
        try
        try 
        cd(dirFly(end).name);
        catch
            try
            cd('extracted');
            dirFly = dir('*fly*extraction');
            cd(dirFly(end).name);
            catch
                dirExtraction = dir('*extraction');
                cd(dirExtraction(1).name);
            end
        end
        dirAnal = dir('analysis*');
        try
           maxFig = openfig([dirAnal(end).folder filesep dirAnal(end).name filesep 'ROI' filesep dirFly(end).name(1:17) 'maxProject.fig']);
        catch
            try
            maxFig = openfig([dirAnal(end).folder filesep dirAnal(end).name filesep 'ROI' filesep dirFly(end).name(1:18) 'maxProject.fig']);
        catch
        maxFig = openfig([dirAnal(end).folder filesep dirAnal(end).name filesep 'ROI' filesep dirFly(end).name(1:19) 'maxProject.fig']);
            end
        end
        drawnow;
        ax1 = gca;
        set(0,'CurrentFigure',plotAllROI);
        s1 = subplot(ceil(length(dirBat)/7),7,day_i);
        fig1 = get(ax1(1),'children');
        copyobj(fig1,s1);
        hold on;
        colormap('gray');
        xlim([0 640]);
        ylim([0 480]);
        set(gca,'CLim',[0 max(max(fig1.CData))/2],'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
        %title(dirBat(day_i).name);
        drawnow;
        cd(dirTop(bat_i).folder);
        close Figure 2;
        catch
            cd(dirTop(bat_i).folder);
        end
    end
    if saveFlag == 1
        saveas(plotAllROI,[saveDir dirBat(bat_i).name(1:2) '_maxProjAllROI_allDays_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotAllROI,[saveDir dirBat(bat_i).name(1:2) '_maxProjAllROI_allDays_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
    end
end
