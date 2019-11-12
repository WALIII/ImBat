function [snakeTracePrePost] = ImBat_plotSnake(snakeTracePrePost,varargin)
% User inputs overrides
saveFlag = 0;

nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'saveflag'
            saveFlag = varargin{i+1};
    end
end


%plot all the flight trajectories concatenated sorted by their preference
snakePlot_prefAll = figure();
imagesc(snakeTracePrePost.normMeanTraceSortPreFlightPost);
colormap(hot);
ylabel('cell number');
set(gca,'yticklabel',{[]});
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
xlabel('time (s)');
sgtitle(['Spatial selectivity sort by flight preference (pre/post/flight): ' snakeTracePrePost.batName ' ' snakeTracePrePost.dateSesh ' ' snakeTracePrePost.sessionType]);
