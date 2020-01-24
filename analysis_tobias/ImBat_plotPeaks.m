function [peaks] = ImBat_plotPeaks(snakeTrace,goodCellIdx,varargin)

batName = snakeTrace.batName;
dateSesh = snakeTrace.dateSesh;
sessionType = snakeTrace.sessionType;
loadFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
saveFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
meanSmooth = 30; %number of video frames to smooth the calcium mean spiking activity

traceDataAll = snakeTrace.meanTracePreFlightPost;
traceDataPre = snakeTrace.meanTracePre;
traceDataFlight = snakeTrace.meanTraceFlight;
traceDataPost = snakeTrace.meanTracePost;

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'batname'
            batName=varargin{i+1};
        case 'datesesh'
            dateSesh = varargin{i+1};
        case 'sessiontype'
            sessionType = varargin{i+1};
        case 'loadflag'
            loadFlag = varargin{i+1};
        case 'saveflag'
            saveFlag = varargin{i+1};
    end
end

%pks = cell(1,length(goodCellIdx));
%locs = cell(1,length(goodCellIdx.goodCellIndex));
%w = cell(1,length(goodCellIdx.goodCellIndex));
%p = cell(1,length(goodCellIdx.goodCellIndex));

smoothAtoB = [];
figure();
for cell_i = goodCellIdx
    smoothAtoB(cell_i,:) = smooth(traceDataAll{1}(cell_i,:),meanSmooth); %smooth the data
    [pks,locs,w,p] = findpeaks(smoothAtoB(cell_i,5:end-5),'MinPeakProminence',0.075);
    plot(traceDataAll{1}(cell_i,5:end-5));
    hold on;
    plot(smoothAtoB(cell_i,5:end-5));
    text(locs+.02,pks,num2str((1:numel(pks))'));
    for pk_i = 1:length(pks)
        tStartAtoMax(cell_i) = locs(pk_i) - length(traceDataPre{1}(pk_i,:)); %find offset from start of flight A to peak
        halfMaxAtoB(pk_i) = pks(pk_i)/2; %find halfmax
        try
            indRiseAtoB(pk_i) = find(smoothAtoB(pk_i,5:locs(pk_i))<=halfMaxAtoB(pk_i),1,'last') + 5; %find where data goes above half max (-5 because we start with the 6th element)
        catch
            indRiseAtoB(pk_i) = 1; %in case the activity never drops below the fwhm by beginning of sample
        end
        try
            indFallAtoB(pk_i) = find(smoothAtoB(cell_i,locs(pk_i):end-5)<=halfMaxAtoB(pk_i),1,'first')+locs(pk_i)-1; %find where data drops below half max (-1 because you need to look 1 back from the first below half max)
        catch
            indFallAtoB(pk_i) = length(smoothAtoB(pk_i,:)); %in case the activity never drops below the fwhm by end of sample
        end
        if indFallAtoB(pk_i) < indRiseAtoB(pk_i)
            indFallAtoB(pk_i) = length(smoothAtoB(pk_i,:));
        end
        fwhmAtoB(pk_i) = indRiseAtoB(pk_i) - indFallAtoB(pk_i); %fwhm by subtracting the indices
        plot(locs(pk_i),pks(pk_i)+0.1,'o','MarkerFaceColor','r','MarkerSize',10) %plot peak on
        plot(indRiseAtoB(pk_i):indFallAtoB(pk_i),halfMaxAtoB(pk_i),'o') %plot the half width max height
        
     end
    
    
    pause
    clf
end