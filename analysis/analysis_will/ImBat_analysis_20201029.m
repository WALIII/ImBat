%% take the first, and last detected day, take the S matrix of first and last cell  (STD)
function out = ImBat_analysis_20201029(CombinedROI,flightPaths,clst)

peak_thresh = 1.5; % SNR threshold
counter = 1;
% find first day
[FlightAlignedROI] = ImBat_Align_FC(CombinedROI,flightPaths,clst);

% %% ROI Analysis
CutCells = FlightAlignedROI.S;
Pre_bound2use = 120; % 2 seond pad
Post_boud2use = 120;
% get dates:
transition_points = find((diff(FlightAlignedROI.CutCells_date)>=1));  %?
transition_points = [1 transition_points size(CutCells,3)];





% get bound to plot
%bound2plot = 1:500;
ROI_ON = FlightAlignedROI.ROI_ON;

bound2plot = round(FlightAlignedROI.ROI_ON-Pre_bound2use:round(FlightAlignedROI.ROI_ON+round(mean(FlightAlignedROI.FlightLength))/120*30)+Post_boud2use);
%ROI_ON = round(FlightAlignedROI.ROI_ON-bound2use);

counter = 1;
col = hsv(size(transition_points,2)+1);


%% ROI Analysis
CutCells = FlightAlignedROI.C_raw;

figure();
hold on;
counter = 1;
for i = 1:size(CutCells,1)
    % draw date lines:
    for ii = 1:size(transition_points,2)-1
        
        adata = zscore(squeeze(CutCells(i,bound2plot,transition_points(ii):transition_points(ii+1))),[],1)'+counter*4;
        [a b] = max(adata');
        var2save(i,ii) = std(b)/sqrt(length(adata)); % cell, transition)
        mn = nanmean(adata); % cell, transition)
        mn = mn-min(mn);
        [a2 b2] = max(mn');
        meanPeakTime(i,ii) = b2;
        meanPeakPeak(i,ii) = a2;
        
    end
    
end

% figure();
% for i = 1:size(CutCells,1);
% cell2plot = i;
%  x = meanPeak(cell2plot,:);
%     y = 1:size(x,2);
%     err = var2save(cell2plot,:)*2;
% errorbar(x,y,err,'horizontal','o')
% xlim([ 0 400]);
% ylim([0 6]);
% pause();
% end

% get the max
% get index of first and last that are not zero
% idx = find(mean(Mean2save)=| 0);
figure();
hold on;
counter = 1;
col = [1,0,0; 0,0,1];

% create sort:
counter = 1;
for i = 1:size(CutCells,1);
    x = meanPeakTime(i,:);
    idx = find(x>1);
    try
        x2sort(i) = [x(idx(1))];
    catch
        x2sort(i) = [x(1)];
        disp('no tracking on any day...');
        idx2remove(counter) = i;
        counter = counter+1;
    end
    
end

if exist('idx2remove','var');
    disp('removing bad indexes');
    CutCells(idx2remove,:,:) = [];
    x2sort(idx2remove) = [];
    meanPeakTime(idx2remove,:) = [];
    meanPeakPeak(idx2remove,:) = [];
end

[a sort_idx] = sort(x2sort);

counter = 1; % reset counter
for i = 1:size(CutCells,1);
    
    cell2plot = sort_idx(i);
    x = meanPeakTime(cell2plot,:);
    err = var2save(cell2plot,:)*4;
    idx = find(x>1);
    y2use = [i i];
    % here is where you decide what day to compare too:
    
    if size(idx,2)>1; 
    x2use = [x(idx(1)) x(idx(2))]; err2use = [err(idx(1)) err(idx(2))]; else x2use = [x(idx(1)) x(idx(end))]; err2use = [err(idx(1)) err(idx(end))]; end
   % x2use = [x(idx(1)) x(idx(end))]; err2use = [err(idx(1)) err(idx(end))]; end
    
    if meanPeakPeak(cell2plot,idx(1)) <peak_thresh; % in
        continue
    end
    for ii = 1:2
        x1 = x2use(ii);
        y1 = y2use(ii);
        err1 = err2use(ii);
        errorbar(x1,y1,err1,'horizontal','o','color',col(ii,:),'CapSize',1);
        % save data
        x_2save(counter,ii) = x2use(ii);
        y_2save(counter,ii) = y2use(ii);
        err_2save(counter,ii) = err2use(ii);
        
    end
    counter = counter+1;
    
end
xlim([ 0 bound2plot(end)]);
ylim([0 i+1]);

% Set where ticks will be
% Get axis handle
ax = gca;
ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
% Set TickLabels;
ylabel('filghts')
xlabel('time from takeoff');
ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};

% plot flight path on top:
ff = FlightAlignedROI.FlightLength;
fl_length =  mean(ff)/120;
fl_length = round(fl_length*30);

plot([ROI_ON ROI_ON],[0 i+1],'--k')
plot([ROI_ON+fl_length ROI_ON+fl_length],[0 i+1],'--k')


% Export data
try
out.fl_length = fl_length;
out.ROI_ON = ROI_ON;
out.x_2save = x_2save;
out.y_2save = y_2save;
out.err_2save = err_2save;
out.fl_length = fl_length;
catch
    disp(' no data exported at this thireshold value')
    out.fl_length = fl_length;
out.ROI_ON = ROI_ON;
out.x_2save = [];
out.y_2save = [];
out.err_2save = [];
end
    
