function out = ImBat_PlotAlignedROIs(FlightAlignedROI,varargin);
% ImBat_PlotAlignedROIs

% Plot individual ROIs aligned to flights from FlightAlignedROI, Rde lines
% indicate dates

% To generate aligned ROIs, first run:
%[FlightAlignedROI] = ImBat_Align_FC(CombinedROI,flightPaths,clst);


% WAL3
% 01/10/2021
display_wait = 1; % pause and display
nparams=length(varargin);
manIn = 0; % manual inputs
manSort = 0; % manual inputs
smooth_data = 0;
% User input
for i=1:2:nparams
    switch lower(varargin{i})
        case 'cells2use' % downsample video
            manIn =1;
            cells2use=varargin{i+1};
        case 'deinterlace' % de-interlace 60fps video
            deinterlace =varargin{i+1};
        case 'sort' % also sort based on input ( like light/dark flights)
            manSort = 1;
            input_order=varargin{i+1};
        case 'smooth' % also sort based on input ( like light/dark flights)
            smooth_data = 1;
            smooth_factor=varargin{i+1};
        case 'display' % also sort based on input ( like light/dark flights)
            display_wait=varargin{i+1};
            
    end
end


figure();
hold on;


% %% ROI Analysis
CutCells = FlightAlignedROI.C;
ROI_ON = FlightAlignedROI.ROI_ON;

% smooth data
if smooth_data ==1;
    disp(['Smoothing data by a factor of', num2str(smooth_factor)]);
    for i = 1:size(CutCells,1)
        for ii = 1:size(CutCells,3)
            CutCells(i,:,ii) = medfilt1(CutCells(i,:,ii),smooth_factor);
        end
    end
end


% get dates:
transition_points = find((diff(FlightAlignedROI.CutCells_date)>=1));  %?
transition_points = [1 transition_points size(CutCells,3)];


% get bound to plot
bound2plot = 1:500;
counter = 1;
col = hsv(size(transition_points,2)+1);
figure();
if manIn ==1;
else
    cells2use = 1:size(CutCells,1);
end
for i = 1:length(cells2use)
    i = cells2use(i);
    subplot(10,1,1:2);% plot shadding for each cell, over lay for each day
    hold on;
    for ii = 1:size(transition_points,2)-1
        
        adata = zscore(squeeze(CutCells(i,bound2plot,transition_points(ii):transition_points(ii+1))),[],1)';
        
        L = size(adata,2);
        se = std(adata)/2;%/10;%sqrt(length(adata));
        mn = nanmean(adata);
        h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col(ii,:)); alpha(0.5);
        plot(mn,'Color',col(ii,:));
    end
    hold off
    subplot(10,1,3:10);
    
    adata = zscore(squeeze(CutCells(i,bound2plot,:)),[],1)';
    %     adata = adata-median(adata(:,500:600)')';
    imagesc(adata,[-3 4]);
    % draw date lines:
    colormap(parula.^2)
    hold on;
    for ii = 1:size(transition_points,2)
        line([1,size(CutCells,2)], [transition_points(ii),transition_points(ii)], 'Color','r','LineWidth',2);
    end
    % Get axis handle
    ax = gca;
    % Set where ticks will be
    ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
    % Set TickLabels;
    ylabel('filghts')
    xlabel('time from takeoff');
    ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};
    title([ 'ROI ' num2str(i)]);
    if display_wait ==1
        pause();
    end
    clf('reset');
    clear adata
end


if manSort ==1; % sort based on manual inputs:
    counter = 1;
    CutCells_1 = CutCells(:,:,input_order);
    CutCells_2 = CutCells;
    CutCells_2(:,:,input_order) = [];
    clear CutCells transition_points
    CutCells = cat(3,CutCells_1,CutCells_2);
    
    transition_points = length(input_order);
    transition_points = [1 transition_points size(CutCells,3)];
    % get bound to plot
    bound2plot = 1:500;
    counter = 1;
    col = hsv(size(transition_points,2)+1);
    figure();
    if manIn ==1;
    else
        cells2use = 1:size(CutCells,1);
    end
    for i = 1:length(cells2use)
        i = cells2use(i);
        subplot(10,1,1:2);% plot shadding for each cell, over lay for each day
        hold on;
        for ii = 1:size(transition_points,2)-1
            
            adata = zscore(squeeze(CutCells(i,bound2plot,transition_points(ii):transition_points(ii+1))),[],1)';
            
            L = size(adata,2);
            se = std(adata)/4;%/10;%sqrt(length(adata));
            mn = nanmean(adata);
            h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col(ii,:)); alpha(0.5);
            plot(mn,'Color',col(ii,:));
        end
        hold off
        subplot(10,1,3:10);
        
        adata = zscore(squeeze(CutCells(i,bound2plot,:)),[],1)';
        out.ROI_sorted{counter} = adata;
        out.transition = transition_points;
        counter = counter+1;
        %     adata = adata-median(adata(:,500:600)')';
        imagesc(adata,[-3 4]);
        % draw date lines:
        colormap(parula.^2)
        hold on;
        for ii = 1:size(transition_points,2)
            line([1,size(CutCells,2)], [transition_points(ii),transition_points(ii)], 'Color','r','LineWidth',2);
        end
        % Get axis handle
        ax = gca;
        % Set where ticks will be
        ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
        % Set TickLabels;
        ylabel('filghts')
        xlabel('time from takeoff');
        ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};
        title([ 'ROI ' num2str(i)]);
        if display_wait ==1
            pause();
        end
        clf('reset');
        clear adata
    end
end






