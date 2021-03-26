function  ScoreMatrix =ImBat_BuildScoreMatrix(CombinedROI,FlightAlignedROI)

% ROI and flight stats
%clst = 2% built from: ImBat_analysis_10212026(flightPaths,ROI_Data,CombinedROI,clst);

%[FlightAlignedROI] = ImBat_Align_FC(CombinedROI,flightPaths,clst);
plot_things = 0;

stop_time = 1;
sub2plot = 5;
TrackScoreThresh = 0.75;

disp ('performing pairwise comparisons');


% %% ROI Analysis
CutCells = FlightAlignedROI.C;
ROI_ON = FlightAlignedROI.ROI_ON;

% get dates:
transition_points = find((diff(FlightAlignedROI.CutCells_date)>=1));  %?
transition_points = [1 transition_points size(CutCells,3)];

% %% ROI Analysis
max_date = max(FlightAlignedROI.CutCells_date);

% get bound to plot
bound2plot = round(FlightAlignedROI.ROI_ON-30:FlightAlignedROI.ROI_ON+mean(FlightAlignedROI.FlightLength)/120*30+60);
trial_cutoff = 10; % min trial number accepted;
counter = 1;
col = hsv(size(transition_points,2)+1);
counter = 1
figure(); 
for i = 1:max_date; % for every day, compair pairs
    for i2 = 1:max_date % pairwise comparison.. 
    % get index for first day:
   idx2use_d1 =  find(FlightAlignedROI.CutCells_date==i);
   idx2use_d2 =  find(FlightAlignedROI.CutCells_date==i2);

   for ii = 1:size(CutCells,1); %for every cell
       % first, check tracking score
       try
       Trackscore = CombinedROI.p_same_registered_pairs{ii}(i,i2);
       catch
           Trackscore = 1; % if we sim data this will exceede bounds...
       end
if  isnan(Trackscore) || Trackscore<TrackScoreThresh
    ScoreMatrix(1,counter) = NaN;
    ScoreMatrix(2,counter) = NaN;
    ScoreMatrix(3,counter) = NaN;
    ScoreMatrix(4,counter) = NaN;
    ScoreMatrix(5,counter) = NaN;
else

adata_d1 = zscore(squeeze(CutCells(ii,bound2plot,idx2use_d1)),[],1)';
adata_d2 = zscore(squeeze(CutCells(ii,bound2plot,idx2use_d2)),[],1)';

%TO DO add track score cutoff
if size(idx2use_d1,2)<trial_cutoff || size(idx2use_d2,2)< trial_cutoff;
 ScoreMatrix(1,counter) = NaN;
 ScoreMatrix(2,counter) = NaN;
 ScoreMatrix(3,counter) = NaN;
 ScoreMatrix(4,counter) = NaN;
 ScoreMatrix(5,counter) = NaN;
else    
M_adata_d1 = smooth(mean(adata_d1,1),30)'; % day 1 mean
M_adata_d2 = smooth(mean(adata_d2,1),30)'; % day 2 mean
    
% within day corr:
M_adata_d1e = mean(adata_d1(1:2:end,:),1); % day 1 even mean
M_adata_d1o = mean(adata_d1(2:2:end,:),1); % day 1 odd mean
M_adata_d2e = mean(adata_d2(1:2:end,:),1); % day 2 even  mean
M_adata_d2o = mean(adata_d2(2:2:end,:),1); % day 2 odd mean

WiDC_d1 = corr(M_adata_d1e',M_adata_d1o'); % compare even odd, day 1
WiDC_d2 = corr(M_adata_d2e',M_adata_d2o'); % compare even odd, day 2;

[peak_mag_adata_d1, peak_M_adata_d1] = max(M_adata_d1);
[min_dat, peak_M_adata_d1] = min(M_adata_d1);
peak_mag_adata_d1 = peak_mag_adata_d1-min_dat;
[scrap peak_M_adata_d2] = max(M_adata_d2);

try
AcDC = corr(M_adata_d1',M_adata_d2'); % compare mean across days
catch
    disp('whoops');
end
ScoreMatrix(1,counter) = Trackscore; % spatial footprint tracking score
ScoreMatrix(2,counter) = AcDC; % across day correlation 
ScoreMatrix(3,counter) = WiDC_d1; % within day correlation, d1 
ScoreMatrix(4,counter) = WiDC_d2; % within day correlation, d2 
ScoreMatrix(5,counter) = abs(i-i2); % seperation across days;
ScoreMatrix(6,counter) = peak_M_adata_d1; % peak, in time of transient, day1
ScoreMatrix(7,counter) = peak_M_adata_d2; % peak, in time of transient, day2
ScoreMatrix(8,counter) = peak_mag_adata_d1; % peak magnitude, day1


counter = counter+1;
end

   end
   end
    end
end


if plot_things ==1;
figure(); 
hold on;
plot(ScoreMatrix(1,:),ScoreMatrix(2,:),'*');
plot([0 1],[0 1],'.-r')
   
% for every transition

% for every cell

% 1. check if cell is detected
% 2. check if cell is detected the next day
% 3. internal check, across day check, save this
plotting = plot_things;

if plotting == 1
% do some quantitative analysis
for i = 1: size(CombinedROI.p_same_registered_pairs,2); GG(:,:,i) = CombinedROI.p_same_registered_pairs{i}; end
GG1 = GG(:);
GG1(isnan(GG1))=[];
figure(); histogram(GG1); title('transition tracking quality'); xlabel('quality'); ylabel('count');


for i = 1: size(CombinedROI.p_same_registered_pairs,2); GG(:,:,i) = CombinedROI.p_same_registered_pairs{i}; end





figure(); hold on;  plot(ScoreMatrix(1,:),ScoreMatrix(2,:),'*'); plot([0 1],[0 1],'.-r')

ix1 = find(ScoreMatrix(5,:)==1);
ix2 = find(ScoreMatrix(5,:)==2);
ix3 = find(ScoreMatrix(5,:)==4);

figure(); 
subplot(1,2,1);

hold on;
plot(sort(ScoreMatrix(2,ix1(1:size(ix3,2)))),'r');
plot(sort(ScoreMatrix(2,ix2(1:size(ix3,2)))),'g');
plot(sort(ScoreMatrix(2,ix3(1:size(ix3,2)))),'b');
title('CDF of across-day temporal corr. for putative pairs ( r = 1 g = 3 b = 5)');
xlabel('ROI pairs');
ylabel('count ( subsampled from max of day 5)')

subplot(1,2,2);
hold on;
plot(sort(ScoreMatrix(1,ix1(1:size(ix3,2)))),'r');
plot(sort(ScoreMatrix(1,ix2(1:size(ix3,2)))),'g');
plot(sort(ScoreMatrix(1,ix3(1:size(ix3,2)))),'b');
title('CDF of across-day spatial corr. for putative pairs ( r = 1 g = 3 b = 5)');
xlabel('ROI pairs');



% withinday vs across day ratio: ( closer to one is good.... means no
change_vector = 1-abs( ((ScoreMatrix(3,:)+ScoreMatrix(4,:))/2) ./ScoreMatrix(2,:));
figure()
hold on;
subplot(1,2,1);
hold on;
plot(sort(change_vector(ix1(1:size(ix3,2)))),'r');
plot(sort(change_vector(ix2(1:size(ix3,2)))),'g');
plot(sort(change_vector(ix3(1:size(ix3,2)))),'b');
title('CDF of ROI stability: [ 1 - abs(within/across) ]');
xlabel('ROI pairs');
ylabel('Stability Score ( 1 = perfect stability) ')
subplot(1,2,2);
hold on;
histogram(change_vector(ix1),'BinWidth',1,'FaceColor','r','Normalization','probability');
histogram(change_vector(ix2),'BinWidth',1,'FaceColor','g','Normalization','probability');
histogram(change_vector(ix3),'BinWidth',1,'FaceColor','b','Normalization','probability');


% sorted by time-series correlation, are there trends in the spatial

[s idx2] = sort(ScoreMatrix(3,:));
figure(); plot(ScoreMatrix(1,:),'*');

% split the data into the highest, and lowest temporal scores

high_temp = find(ScoreMatrix(1,:)>=0.9);
low_temp = find(ScoreMatrix(1,:)<=0.1);

figure(); 
hold on;
histogram((ScoreMatrix(2,high_temp)),'BinWidth',0.2,'FaceColor','b','Normalization','probability');
histogram((ScoreMatrix(2,low_temp)),'BinWidth',0.2,'FaceColor','r','Normalization','probability');
xlabel('trajectory correlation');
ylabel('frequency');
title(' top 10 (blue) and bottem 10% (red) spatial correlations');




% Plot correlation over days 
clear Var1 Var2
for i = 1:max(ScoreMatrix(5,:));
    
ix1 = find(ScoreMatrix(5,:)==i & ScoreMatrix(1,:)< 0.1);

Var1(i) = nanmean(ScoreMatrix(2,ix1));
Var2(i) = nanstd(ScoreMatrix(2,ix1));
end

% concat within day stability on day 1:
ix1 = find(ScoreMatrix(5,:)==1);
toAdd = ScoreMatrix(3,ix1);
Var1 = cat(2,nanmean(toAdd),Var1);
Var2 = cat(2,nanstd(toAdd),Var2);

figure();
errorbar([1:size(Var1,2)], Var1,...
        (Var2)/sqrt(length(Var2)), 'color', ...
     'b', 'LineWidth', 2);
 title('Tuning stability over days,all ROIs');
 xlabel('days');
 ylabel('tuning stability to day 1');

end
 

plotting2 = 0;
if plotting2 ==1;
% get bound to plot
bound2plot = 1:400;
counter = 1;
col = hsv(size(transition_points,2)+1);
figure();
for i = 1:size(CutCells,1);
subplot(1,2,1);
adata = zscore(squeeze(CutCells(i,bound2plot,:)),[],1)';
    imagesc(adata);
    % draw date lines:
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
    
    subplot(1,2,2);
    h = CombinedROI.p_same_registered_pairs{i};
     h2 = nanmean(h);
     idx2use = 1-isnan(h2);
    h(isnan(h)) = 1;
    h = h.*idx2use;
    h = h.*idx2use';
imagesc(h,[-0.1 1]); colorbar
    pause();
    clf('reset');
    clear adata
end
end


end

