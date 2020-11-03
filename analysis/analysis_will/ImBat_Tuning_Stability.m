
function ImBat_Tuning_Stability(ScoreMatrix);
% Plot data stability over days:

% run:
% ScoreMatrix = ImBat_analysis_10212026(flightPaths,ROI_Data,CombinedROI,3);

% Plot correlation over days 
clear Var1 Var2
for i = 1:max(ScoreMatrix(5,:));
    
ix1 = find(ScoreMatrix(5,:)==i);

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

 
 
 % top and bottem confidence scores
 
% Plot correlation over days 
clear Var1 Var2 Var3 Var4
for i = 1:max(ScoreMatrix(5,:));
    
ix1 = find(ScoreMatrix(5,:)==i & ScoreMatrix(1,:)> 0.9);
ix2 = find(ScoreMatrix(5,:)==i & ScoreMatrix(1,:)< 0.2);

Var1(i) = nanmean(ScoreMatrix(2,ix1));
Var2(i) = nanstd(ScoreMatrix(2,ix1));

Var3(i) = nanmean(ScoreMatrix(2,ix2));
Var4(i) = nanstd(ScoreMatrix(2,ix2));
end

% concat within day stability on day 1:
ix1 = find(ScoreMatrix(5,:)==i & ScoreMatrix(1,:)> 0.95);
ix2 = find(ScoreMatrix(5,:)==i & ScoreMatrix(1,:)< 0.2);

toAdd = ScoreMatrix(3,ix1);
toAdd3 = ScoreMatrix(3,ix2);

Var1 = cat(2,nanmean(toAdd),Var1);
Var2 = cat(2,nanstd(toAdd),Var2);
Var3 = cat(2,nanmean(toAdd3),Var3);
Var4 = cat(2,nanstd(toAdd3),Var4);


figure();
hold on;
errorbar([1:size(Var1,2)], Var1,...
        (Var2)/sqrt(length(Var2)), 'color', ...
     'r', 'LineWidth', 2);
 
 errorbar([1:size(Var3,2)], Var3,...
        (Var4)/sqrt(length(Var4)), 'color', ...
     'g', 'LineWidth', 2);
 title('Tuning stability over days, for best (r) and worst ( g) tracked cells');
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

 
 %% take the first, and last detected day, take the S matrix of first and last cell  (STD)
 
 % find 
 
 
 %% Stability 
 
 