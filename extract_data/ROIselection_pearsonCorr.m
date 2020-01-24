function [cellsPearsonCorr,goodCellIdx] = ROIselection_pearsonCorr(snakeTrace,goodCellIdx,varargin)

traceData = snakeTrace.tracePreFlightPost; %which data stream to use?
batName = snakeTrace.batName;
dateSesh = snakeTrace.dateSesh;
sessionType = snakeTrace.sessionType;
saveFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
highCorrThresh = 0.5; %threshold for which cells are deemed to be high correlation to use for analyses
plotOddEvenAllFlag = 1;
plotHighCorrFlag = 1;
plotHeatMapFlag = 1;

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'saveflag'
            saveFlag = varargin{i+1};
    end
end
%% calculate the mean neural traces for odd and even flights (& normalize)

%initialize the variables for odd and even flights
meanTraceFlightOdd = zeros(length(traceData{2}(1,1,:)),length(traceData{2}(1,:,1)));
meanTraceFlightEven = zeros(length(traceData{2}(1,1,:)),length(traceData{2}(1,:,1)));
normMeanTraceFlightOdd = zeros(length(traceData{2}(1,1,:)),length(traceData{2}(1,:,1)));
normMeanTraceFlightEven = zeros(length(traceData{2}(1,1,:)),length(traceData{2}(1,:,1)));
traceFlightOdd = zeros(floor(length(traceData{2}(:,1,1))/2),length(traceData{2}(1,:,1)),length(traceData{2}(1,1,:)));
traceFlightEven = zeros(floor(length(traceData{2}(:,1,1))/2),length(traceData{2}(1,:,1)),length(traceData{2}(1,1,:)));
for cell_i = goodCellIdx.goodCellIndex%length(traceData{2}(1,1,:))
    %build vector of odd flights for each cell
    for flight_i = 1:2:length(traceData{2}(:,1,1))
        traceFlightOdd(flight_i,:,cell_i) = traceData{2}(flight_i,:,cell_i);
        %traceFlightOdd(flight_i,:,cell_i) = traceFlightOdd(flight_i,:,cell_i) - abs(min(traceData{2}(flight_i,:,cell_i)));
    end
    %build vector of even flights for each cell
    for flight_i = 2:2:length(traceData{2}(:,1,1))
        traceFlightEven(flight_i,:,cell_i) = traceData{2}(flight_i,:,cell_i);
        %traceFlightEven(flight_i,:,cell_i) = traceFlightEven(flight_i,:,cell_i) - abs(min(traceData{2}(flight_i,:,cell_i)));
    end
    %calculate the means for all the odd/even flights and normalize these traces
    meanTraceFlightOdd(cell_i,:) = mean(traceFlightOdd(1:2:end,:,cell_i),1);
    meanTraceFlightEven(cell_i,:) = mean(traceFlightEven(2:2:end,:,cell_i),1);
    normMeanTraceFlightOdd(cell_i,:) = zscore(meanTraceFlightOdd(cell_i,:));
    normMeanTraceFlightOdd(cell_i,:) = normMeanTraceFlightOdd(cell_i,:) - min(normMeanTraceFlightOdd(cell_i,:));
    normMeanTraceFlightEven(cell_i,:) = zscore(meanTraceFlightEven(cell_i,:));
    normMeanTraceFlightEven(cell_i,:) = normMeanTraceFlightEven(cell_i,:) - min(normMeanTraceFlightEven(cell_i,:));
end
%perform pearson correlation of the odd vs even flights
[rho,pval] = corr(meanTraceFlightOdd',meanTraceFlightEven');

%plot the histogram of all cells vs all cells correlation
histPearsonCorr = figure();
histogram(rho);
p1=subplot(2,1,1);
histogram(rho);
title('Pearson Corr: odd vs even trials, all neurons (r)');
p2=subplot(2,1,2);
histogram(pval);
title('Pearson Corr: odd vs even trials, all neurons (pval)');
sgtitle([snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);

%build the vector of the r and pval for only the cell-to-cell correlations
for cell_i = 1:length(rho(1,:))
    pairRho(cell_i) = rho(cell_i,cell_i);
    pairPval(cell_i) = pval(cell_i,cell_i);
end
%plot the histogram of paired correlation values
histPairedCorr = figure();
p1=subplot(2,1,1);
histogram(pairRho,10);
title('Pearson Corr Odd vs Even Trials (r)');
p2=subplot(2,1,2);
histogram(pairPval);
title('Pearson Corr Odd vs Even Trials (pval)');
sgtitle([snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);

if saveFlag ==1
    if ~exist([pwd,'/PearsonCorr'])>0;
        
        mkdir('PearsonCorr');
    end
    % Save 'each cell' as jpg and fig files..
    set(findall(histPairedCorr,'-property','FontSize'),'FontSize',20);
    saveas(gcf,[pwd '/PearsonCorr/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_hist_PearsonCorr_pairedCells30sec.tif']);
    savefig(gcf,[pwd '/PearsonCorr/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_hist_PearsonCorr_pairedCells30sec.fig']);
    saveas(gcf, [pwd '/PearsonCorr/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_hist_PearsonCorr_pairedCells30sec.svg']);
    set(findall(histPearsonCorr,'-property','FontSize'),'FontSize',20);
    saveas(gcf,[pwd '/PearsonCorr/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_hist_PearsonCorr_pairedCells30sec.tif']);
    savefig(gcf,[pwd '/PearsonCorr/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_hist_PearsonCorr_pairedCells30sec.fig']);
    saveas(gcf, [pwd '/PearsonCorr/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_hist_PearsonCorr_pairedCells30sec.svg']);
    
    
end
%% plot the means and norm means for all the cells along with the pearson corrs
if plotOddEvenAllFlag ==1
    plotOddEvenFlights = figure('units','normalized','outerposition',[0 0 1 0.5]);
    for cell_i = goodCellIdx.goodCellIndex%length(traceData{2}(1,1,:))
        %plot all of the odd traces in left figure
        p1 = subplot(1,3,1);
        for flight_i = 1:2:length(traceData{2}(:,1,1))
            plot(traceFlightOdd(flight_i,:,cell_i));
            hold on;
        end
        title(['Odd flights: ROI# ' num2str(cell_i)]);
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        ylabel('df/f');
        %plot all of the even traces in right figure
        p2 = subplot(1,3,2);
        for flight_i = 2:2:length(traceData{2}(:,1,1))
            plot(traceFlightEven(flight_i,:,cell_i));
            hold on;
        end
        title(['Even flights: ROI# ' num2str(cell_i)]);
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        ylabel('df/f');
        %plot the mean traces for odd and even in same plot +normalized
        p3 = subplot(1,3,3);
        plot(meanTraceFlightOdd(cell_i,:),'Color','r');
        hold on
        plot(meanTraceFlightEven(cell_i,:),'Color','b');
        plot(normMeanTraceFlightOdd(cell_i,:),'Color','g');
        plot(normMeanTraceFlightEven(cell_i,:),'Color','k');
        legend('Odd','Even','Norm Odd','Norm Even');
        title(['Means of odd and even flights: ROI# ' num2str(cell_i)]);
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(1)+0.05*xlim(2),ylim(1)+0.07*ylim(2),{['r= ' num2str(pairRho(cell_i))];['pval =' num2str(pairPval(cell_i))]});
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        ylabel('df/f');
        sgtitle([snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
        
        if saveFlag ==1
            if ~exist([pwd,'/PearsonCorr'])>0
                
                mkdir('PearsonCorr');
            end
            % Save 'each cell' as jpg and fig files..
            %set(findall(gcf,'-property','FontSize'),'FontSize',20);
            saveas(gcf,[pwd '/PearsonCorr/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_meanTrace_oddEven_PearsonCorr30sec-' num2str(cell_i) '.tif']);
            savefig(gcf,[pwd '/PearsonCorr/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_meanTrace_oddEven_PearsonCorr30sec-' num2str(cell_i) '.fig']);
            saveas(gcf, [pwd '/PearsonCorr/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_meanTrace_oddEven_PearsonCorr30sec-' num2str(cell_i) '.svg']);
        else
            pause
        end
        clf;
    end
    close(plotOddEvenFlights)
end
%% plot the mean and normalized odd/even traces for cells above high correlation threshold
if plotHighCorrFlag ==1
    %initialize the variables to store the indices of high and low correlation cells
    highCorrIndex = [];
    lowCorrIndex = [];
    
    plotHighCorr = figure();
    for cell_i = 1:length(rho(1,:))
        if pairRho(cell_i)>=highCorrThresh
            highCorrIndex = [highCorrIndex cell_i];
            plot(meanTraceFlightOdd(cell_i,:),'Color','r');
            hold on
            plot(meanTraceFlightEven(cell_i,:),'Color','b');
            plot(normMeanTraceFlightOdd(cell_i,:),'Color','g');
            plot(normMeanTraceFlightEven(cell_i,:),'Color','k');
            legend('Odd','Even','Norm Odd','Norm Even');
            title([snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType ': High Corr odd and even flights: ROI# ' num2str(cell_i)]);
            ylim=get(gca,'ylim');
            xlim=get(gca,'xlim')
            text(xlim(1)+0.05*xlim(2),ylim(1)+0.07*ylim(2),{['r= ' num2str(pairRho(cell_i))];['pval =' num2str(pairPval(cell_i))]});
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
            ylabel('df/f');
            
            if saveFlag ==1
                if ~exist([pwd,'/PearsonCorr/highCorr_cells'])>0
                    mkdir([pwd,'/PearsonCorr/highCorr_cells']);
                end
                % Save 'each cell' as jpg and fig files..
                %set(findall(gcf,'-property','FontSize'),'FontSize',20);
                saveas(gcf,[pwd '/PearsonCorr/highCorr_cells/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_meanTrace_oddEven_highCorr30sec-' num2str(cell_i) '.tif']);
                savefig(gcf,[pwd '/PearsonCorr/highCorr_cells/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_meanTrace_oddEven_highCorr30sec-' num2str(cell_i) '.fig']);
                saveas(gcf, [pwd '/PearsonCorr/highCorr_cells/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_meanTrace_oddEven_highCorr30sec-' num2str(cell_i) '.svg']);
            else
                pause
            end
            clf
        elseif pairRho(cell_i)<highCorrThresh
            lowCorrIndex = [lowCorrIndex cell_i];
        end
    end
    close(plotHighCorr)
end
%% Save the structure
cellsPearsonCorr.rho = rho;
cellsPearsonCorr.pval = pval;
cellsPearsonCorr.pairRho = pairRho;
cellsPearsonCorr.pairPval = pairPval;
cellsPearsonCorr.meanTraceFlightEven = meanTraceFlightEven;
cellsPearsonCorr.meanTraceFlightOdd = meanTraceFlightOdd;
cellsPearsonCorr.normMeanTraceFlightEven = normMeanTraceFlightEven;
cellsPearsonCorr.normMeanTraceFlightOdd = normMeanTraceFlightOdd;
cellsPearsonCorr.traceFlightEven = traceFlightEven;
cellsPearsonCorr.traceFlightOdd = traceFlightOdd;
cellsPearsonCorr.highCorrThresh = highCorrThresh;
if plotHighCorrFlag == 1
cellsPearsonCorr.highCorrIndex = highCorrIndex;
goodCellIdx.highCorrIndex = highCorrIndex;
end
goodCellIdx.highCorrThresh = highCorrThresh;
if saveFlag ==1
    save([pwd '/' batName '_' dateSesh '_' sessionType '_cellsPearsonCorr30sec.mat'],'cellsPearsonCorr');
    save([pwd '/' batName '_' dateSesh '_' sessionType '_goodCellIdx.mat'],'goodCellIdx');
end

%% plot the matrix of odd vs even flights (sorted by the odd flights)
if plotHeatMapFlag ==1
[~,maxOdd] = max(cellsPearsonCorr.meanTraceFlightOdd,[],2); 
[Bodd,Iodd] = sort(maxOdd);
oddSortedHeatMap = figure(); 
subplot(2,1,1); 
imagesc(cellsPearsonCorr.meanTraceFlightOdd(Iodd,:)); 
title('Odd'); 
subplot(2,1,2); 
imagesc(cellsPearsonCorr.meanTraceFlightEven(Iodd,:));%,[10 50]); 
title('Even'); 
colormap('hot');
sgtitle(['Odd and even flights sorted by odd: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);

if saveFlag ==1
    if ~exist([pwd,'/PearsonCorr'])>0
        mkdir('PearsonCorr');
    end
    % Save 'each cell' as jpg and fig files..
    %set(findall(gcf,'-property','FontSize'),'FontSize',20);
    saveas(oddSortedHeatMap,[pwd '/PearsonCorr/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_oddSortedHeatMap30sec.tif']);
    savefig(oddSortedHeatMap,[pwd '/PearsonCorr/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_oddSortedHeatMap30sec.fig']);
    saveas(oddSortedHeatMap, [pwd '/PearsonCorr/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_oddSortedHeatMap30sec.svg']);
end
end
