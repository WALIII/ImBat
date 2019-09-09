date = 'zu190809';
label = 'Zuzu_190809_fly-1';

saveas(flightPathsClusterAll, ['X:\users\tobias\flight\processed\analysis_done_190905\' date '\' label '_extraction\analysis\flights\' label '_flightPathsClusterAll.tif']);
saveas(flightPathsClusterEach, ['X:\users\tobias\flight\processed\analysis_done_190905\' date '\' label '_extraction\analysis\flights\' label '_flightPathsClusterEach.tif']);
save(['X:\users\tobias\flight\processed\analysis_done_190905\' date '\' label '_extraction\analysis\' label '_flightPaths.mat'],'flightPaths')
savefig(flightPathsClusterAll,['X:\users\tobias\flight\processed\analysis_done_190905\' date '\' label '_extraction\analysis\flights\' label '_flightPathsClusterAll.fig']);
savefig(flightPathsClusterEach,['X:\users\tobias\flight\processed\analysis_done_190905\' date '\' label '_extraction\analysis\flights\' label '_flightPathsClusterEach.fig']);



saveas(flightPathsClusterToFeederAll, ['X:\users\tobias\flight\processed\analysis_done_190905\' date '\' label '_extraction\analysis\flights\' label '_flightPathsClusterToFeederAll.tif']);
saveas(flightPathsClusterToFeederEach, ['X:\users\tobias\flight\processed\analysis_done_190905\' date '\' label '_extraction\analysis\flights\' label '_flightPathsClusterToFeederEach.tif']);
save(['X:\users\tobias\flight\processed\analysis_done_190905\' date '\' label '_extraction\analysis\' label '_flightsToFromFeeders.mat'],'flightFeedersStartStop')
savefig(flightPathsClusterToFeederAll,['X:\users\tobias\flight\processed\analysis_done_190905\' date '\' label '_extraction\analysis\flights\' label '_flightPathsClusterToFeederAll.fig']);
savefig(flightPathsClusterToFeederEach,['X:\users\tobias\flight\processed\analysis_done_190905\' date '\' label '_extraction\analysis\flights\' label '_flightPathsClusterToFeederEach.fig']);