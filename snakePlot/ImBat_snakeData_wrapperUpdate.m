function ImBat_snakeData_wrapperUpdate
plotFlag = 0;
saveFlag = 1;
g = dir('G*');
z = dir('Z*');
dirTop = vertcat(g,z); %find all folders in top quality directory

%for each bat/session
for sesh_i = 1:length(dirTop)
    
    %get meta info for each bat/day
    cd([dirTop(sesh_i).name filesep 'extracted']);
    dirFly = dir('*fly*extraction*');
    batName = dirFly(1).name(1:3);
    dateSesh = dirFly(1).name(5:10);
    sessionType = dirFly(1).name(12:16);
    fileName = [batName '_' dateSesh '_' sessionType];
    %load cellData, flightpaths,alignment files
    cd(dirFly(1).name);
    dirAnal = dir('analysis*');
    dirFP = dir([dirAnal(end).folder filesep dirAnal(end).name filesep '*flightPaths.mat']);
    load([dirFP.folder filesep dirFP.name]);
    dirProc = dir('processed*');
    dirCD = dir([dirProc(end).folder filesep dirProc(end).name filesep 'results.mat']);
    dirA = dir([dirProc(end).folder filesep dirProc(end).name filesep 'Alignment.mat']);
    cellData = load([dirCD.folder filesep dirCD.name]);
    alignment = load([dirA.folder filesep dirA.name]);
    
    [snakeTrace_cRaw,snakeTrace_c,snakeTrace_s] = ImBat_snakeData(cellData,flightPaths,alignment,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
    %add the 3 data sets to snakeTrace
    snakeTrace.cRaw = snakeTrace_cRaw;
    snakeTrace.c = snakeTrace_c;
    snakeTrace.s = snakeTrace_s;
    if plotFlag == 1
        [snakeTrace_plots] = ImBat_plotSnake(snakeTrace_cRaw);
        if saveFlag == 1
            if ~exist([dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots'])
                mkdir([dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots']);
            end
            saveas(snakeTrace_plots.snakePlot_clust, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlots_clust.svg']);
            saveas(snakeTrace_plots.snakePlot_clustOddEven, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlots_clustOddEven.svg']);
            saveas(snakeTrace_plots.snakePlot_clustBy1, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlots_clustBy1.svg']);
            saveas(snakeTrace_plots.snakePlot_prefEachPrePostFlight, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlot_prefEachPrePostFlight.svg']);
            saveas(snakeTrace_plots.snakePlot_clustPrePostFlight, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlot_clustPrePostFlight.svg']);
            saveas(snakeTrace_plots.snakePlot_clust, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlots_clustAll.tif']);
            saveas(snakeTrace_plots.snakePlot_clustOddEven, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlots_clustOddEven.tif']);
            saveas(snakeTrace_plots.snakePlot_clustBy1, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlots_clustBy1.tif']);
            saveas(snakeTrace_plots.snakePlot_prefEachPrePostFlight, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlot_prefEachPrePostFlight.tif']);
            saveas(snakeTrace_plots.snakePlot_clustPrePostFlight, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlot_clustPrePostFlight.tif']);
            saveas(snakeTrace_plots.snakePlot_prefAll, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlot_prefAll.tif']);
            saveas(snakeTrace_plots.snakePlot_prefAll, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlot_prefAll.svg']);
            saveas(snakeTrace_plots.snakePlot_clustBy1PrePostFlight, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlot_clustBy1PrePostFlight.tif']);
            saveas(snakeTrace_plots.snakePlot_clustBy1PrePostFlight, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlot_clustBy1PrePostFlight.svg']);
            savefig(snakeTrace_plots.snakePlot_clust, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlots_clustAll.fig']);
            savefig(snakeTrace_plots.snakePlot_clustOddEven, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlots_clustOddEven.fig']);
            savefig(snakeTrace_plots.snakePlot_clustBy1, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlots_clustBy1.fig']);
            savefig(snakeTrace_plots.snakePlot_clustPrePostFlight, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlot_clustPrePostFlight.fig']);
            savefig(snakeTrace_plots.snakePlot_prefEachPrePostFlight, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlot_prefEachPrePostFlight.fig']);
            savefig(snakeTrace_plots.snakePlot_prefAll, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlots_prefAll.fig']);
            savefig(snakeTrace_plots.snakePlot_clustBy1PrePostFlight, [dirAnal(end).folder filesep dirAnal(end).name filesep 'snakePlots' filesep fileName '_snakePlot_clustBy1PrePostFlight.fig']);
        end
    end
    if saveFlag == 1
        save([dirAnal(end).folder filesep dirAnal(end).name filesep fileName '_snakePlotData.mat'],'snakeTrace');
    end
    cd(dirTop(sesh_i).folder);
    close all;
end