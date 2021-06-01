
function ImBat_Tuning_Stability(ScoreMatrix,col);
% Plot data stability over days:

% run:
% ScoreMatrix = ImBat_analysis_10212026(flightPaths,ROI_Data,CombinedROI,3);

% Plot correlation over days
plot_things = 1;
plot_all_things = 0;
plot_error_bars = 1; % shadded error instead if zero

% Variables:

minROI_TrackingScore = 0.5; % only effects green trace
minROI_SNR = 1.5;
for i = 1:max(ScoreMatrix(5,:));
    
    ix1 = find(ScoreMatrix(5,:)==i);
    
    try
        Var1(i) = nanmedian(ScoreMatrix(2,ix1));
        Var2(i) = nanstd(ScoreMatrix(2,ix1));
    catch
        disp('No data for this epoch and settings');
        return
    end
end

% concat within day stability on day 1:
ix1 = find(ScoreMatrix(5,:)==1);
toAdd = ScoreMatrix(3,ix1);
Var1 = cat(2,nanmedian(toAdd),Var1);
Var2 = cat(2,nanstd(toAdd),Var2);
%
% figure();
% errorbar([1:size(Var1,2)], Var1,...
%         (Var2)/sqrt(length(Var2)), 'color', ...
%      'b', 'LineWidth', 2);
%  title('Tuning stability over days,all ROIs');
%  xlabel('days');
%  ylabel('tuning stability to day 1');



% top and bottem confidence scores

% Plot correlation over days
clear Var1 Var2 Var3 Var4
for i = 1:max(ScoreMatrix(5,:));
    
    ix1 = find(ScoreMatrix(5,:)==i & ScoreMatrix(1,:)> minROI_TrackingScore & ScoreMatrix(8,:)> minROI_SNR); % find good tracking, and high peak value
    ix2 = find(ScoreMatrix(5,:)==i & ScoreMatrix(1,:)< 0.2 & ScoreMatrix(8,:)> minROI_SNR);
    ix3 = find(ScoreMatrix(5,:)==i & ScoreMatrix(8,:)> minROI_SNR);
    
    
    Var1(i) = nanmedian(ScoreMatrix(2,ix1));
    Var2(i) = nanstd(ScoreMatrix(2,ix1));
    
    Var3(i) = nanmedian(ScoreMatrix(2,ix2));
    Var4(i) = nanstd(ScoreMatrix(2,ix2));
    
    Var5(i) = nanmedian(ScoreMatrix(2,ix3));
    Var6(i) = nanstd(ScoreMatrix(2,ix3));
end

% concat within day stability on day 1:
ix1 = find(ScoreMatrix(5,:)==i & ScoreMatrix(1,:)> minROI_TrackingScore & ScoreMatrix(8,:)> minROI_SNR); % find good tracking, and high peak value
ix2 = find(ScoreMatrix(5,:)==i & ScoreMatrix(1,:)< 0.2 & ScoreMatrix(8,:)> minROI_SNR);
ix3 = find(ScoreMatrix(5,:)==i & ScoreMatrix(8,:)> minROI_SNR);


toAdd = ScoreMatrix(3,ix1);
toAdd2 = ScoreMatrix(3,ix2);
toAdd3 = ScoreMatrix(3,ix3);

% high corr
Var1 = cat(2,nanmedian(toAdd),Var1);
Var2 = cat(2,nanstd(toAdd),Var2);
% Low Corr
Var3 = cat(2,nanmedian(toAdd2),Var3);
Var4 = cat(2,nanstd(toAdd2),Var4);
% all
Var5 = cat(2,nanmedian(toAdd3),Var5);
Var6 = cat(2,nanstd(toAdd3),Var6);

% Explort data for group plotting
out.Var1 = Var5;
out.Var2 = Var6;


if plot_things ==1;
    
    if plot_error_bars ==1
        hold on;
        errorbar([1:size(Var5,2)], Var5,...
            (Var6)/sqrt(length(Var6)), 'color', ...
            col, 'LineWidth', 2);
        title('Tuning stability over days, for best (r) and worst ( g) and all ( b) tracked cells');
        xlabel('days');
        ylabel('tuning stability to day 1');
        
        % STATS
        % STATS

X1 = [1:size(Var5,2)]';
y = Var5';
x1 = ones(size(Var5,2),1);
X = [x1 X1];    % Includes column of ones
[b,~,~,~,stats] = regress(y,X)


    else
        hold on;
        L = size(Var5,2);
        se = (Var6)/sqrt(length(Var6));%sqrt(length(adata));
        mn = Var5;
        h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col); alpha(0.5);
        plot(mn,'Color',col);
    end
end


if plot_all_things ==1;
    
    % Plot data
    
    figure();
    hold on;
    errorbar([1:size(Var1,2)], Var1,...
        (Var2)/sqrt(length(Var2)), 'color', ...
        'r', 'LineWidth', 2);
    
    errorbar([1:size(Var3,2)], Var3,...
        (Var4)/sqrt(length(Var4)), 'color', ...
        'g', 'LineWidth', 2);
    
    errorbar([1:size(Var5,2)], Var5,...
        (Var6)/sqrt(length(Var6)), 'color', ...
        'b', 'LineWidth', 2);
    title('Tuning stability over days, for best (r) and worst ( g) and all ( b) tracked cells');
    xlabel('days');
    ylabel('tuning stability to day 1');
    
    
    
    
    % plot peaks over time
    
    
    
    
    
    
    
    
    
    
    
    % EXTRA
    
    
    
    
    
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
    
end


%
% % Plot correlation over days
% clear Var1 Var2
% for i = 1:max(ScoreMatrix(5,:));
%
% %ix1 = find(ScoreMatrix(5,:)==i & ScoreMatrix(1,:)> 0.7 & ScoreMatrix(8,:)> 1.2); % find good tracking, and high peak value
% ix1 = find(ScoreMatrix(5,:)==i &  ScoreMatrix(8,:)> 1.2); % find good tracking, and high peak value
%
% Var1(i) = nanmedian(ScoreMatrix(2,ix1));
% Var2(i) = nanstd(ScoreMatrix(2,ix1));
% end
%
% % concat within day stability on day 1:
% ix1 = find(ScoreMatrix(5,:)==1);
% toAdd = ScoreMatrix(3,ix1);
% Var1 = cat(2,nanmedian(toAdd),Var1);
% Var2 = cat(2,nanstd(toAdd),Var2);
%
% figure();
% errorbar([1:size(Var1,2)], Var1,...
%         (Var2)/sqrt(length(Var2)), 'color', ...
%      'b', 'LineWidth', 2);
%  title('Tuning stability over days,all ROIs');
%  xlabel('days');
%  ylabel('tuning stability to day 1');


%% take the first, and last detected day, take the S matrix of first and last cell  (STD)

% find


%% Stability

