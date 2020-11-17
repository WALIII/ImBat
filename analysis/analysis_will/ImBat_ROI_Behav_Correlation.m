function out2 = ImBat_ROI_Behav_Correlation(FlightAlignedROI)

% find if there is a correlation to the ROI and behavioral variance 

counter = 1; 

% NOTE: there may be more than one peak per ROI!

% to get data, run:
for ii = 1:size(FlightAlignedROI{1}.C,1) % first 10 cells:
out =  ImBat_analysis_11112020(FlightAlignedROI,ii);

% plot mean of ROI data
% Zscore data:
AllColZ = zscore(out.data2.AllCol);
%remove zeros
xtemp = mean(AllColZ);
clear h1 h2
h1 = find(xtemp ==0);
h2 = find(xtemp ~=0);
AllColZ(:,h1) = [];
x = mean(AllColZ,2);
% get peaks
[Bpks,Blocs] = findpeaks(x,'MinPeakProminence',.5,'MinPeakDistance',1);

% only use peaks in the flight pad ( exclude reward...
Blocs(Blocs<250) = [];
Blocs(Blocs>1200) = [];

% Blocs(Blocs<300) = [];
% Blocs(Blocs>800) = [];


% figure(); 
% hold on;
% plot(x);
% 
% plot(Blocs,ones(size(Blocs))+1,'*r');

% now, get the distributions for 1. flights 2. calcium;
if isempty(Blocs)
else
for iii = 1: size(Blocs,1); 
calPeaks = out.data2.AllCol(Blocs(iii),:);
flPeaksX = squeeze(out.data2.AllFlight(Blocs(iii),1,:));
flPeaksY = squeeze(out.data2.AllFlight(Blocs(iii),2,:));
flPeaksZ = squeeze(out.data2.AllFlight(Blocs(iii),3,:));
calPeakLoc = Blocs(iii);
% figure(); 
% hold on;
% plot(calPeaks,flPeaksX-mean(flPeaksX),'*r')
% plot(calPeaks,flPeaksY-mean(flPeaksY),'*g')
% plot(calPeaks,flPeaksZ-mean(flPeaksZ),'*b')

% % summ the differce from the mean
% a1 = abs(flPeaksX-mean(flPeaksX))+ abs(flPeaksY-mean(flPeaksY)) + abs(flPeaksZ-mean(flPeaksZ));
% 
% figure(); 
% plot(a1,calPeaks,'*r')

% get  position to mean point 
meanPoint(1) = mean(flPeaksX(h2));
meanPoint(2) = mean(flPeaksY(h2));
meanPoint(3) = mean(flPeaksZ(h2));

for i = 1: size(flPeaksZ,1)
    currPoint(1) = flPeaksX(i);
currPoint(2) = flPeaksY(i);
currPoint(3) = flPeaksZ(i);
D(i) = vecnorm(currPoint - meanPoint, 2, 2);
end

% figure(); histogram(D);


out2.D{counter} = D;
out2.calPeaks{counter} = calPeaks;
out2.calPeakLoc{counter} = calPeakLoc;
counter = counter+1;
end
end

end


%% Plotting: ( will also work outside of function) 

figure();
for i = 1: size(out2.D,2)
    
    % remove zeros
    xx = out2.calPeaks{i};
    h = find(xx ==0);

    F2use = out2.D{i};
    C2use = out2.calPeaks{i};
    F2use(h) = [];
    C2use(h) = [];
    
% normalize calcium from 0-1
    C2use = mat2gray(C2use);
    % ..now take the mean, and variance...
    
    varF(i) = var(F2use);
    varC(i) = std(C2use);
    
    meanF(i) = mean(F2use);
    meanC(i) = mean(C2use);
    
    Loca(i) = out2.calPeakLoc{i};
end
Loca = mat2gray(Loca);

% figure()
% plot(meanF,varC,'*');
% % 
% figure();
% plot(varF,varC,'*');



% ========== stats ( regression) ===========%
clear y x1 X
X1 = meanF';
y = meanC';
x1 = ones(size(X1,1),1);
X = [x1 X1];    % Includes column of ones

[~,~,~,~,stats] = regress(y,X);
X = [X1];
mdl = fitlm(X,y);
figure(); 
plotAdded(mdl)
title('Calcium consistancy vs Flight Variance');
ylabel(' Mean, Normalized Ca+ activity');
xlabel('mean Flight distance (cm);');



% More plotting..... :

figure(); 
plot(meanF,(meanC),'*');
title('Mean Calcium intensity vs aveage flight distance from the mean flight');
ylabel(' Mean, Normalized Ca++ activity');
xlabel(' average euclidian distance from the mean flight in that cluster');



figure();
subplot(2,1,1);

hold on;
plot(Loca,(meanC),'*');
[g1 g2] = sort(Loca);
plot(Loca(g2),(smooth(meanC(g2),30,'lowess')),'LineWidth',3);

title('Mean Calcium intensity vs Flight Phase');
ylabel(' Mean, Normalized Ca++ activity');
xlabel(' Flight Phase');




subplot(2,1,2);
hold on;
plot(Loca,(meanF),'*');
plot(Loca(g2),(smooth(meanF(g2),30,'lowess')),'LineWidth',3);
title('Mean Flight distance from center of cluster vs Flight Phase');
ylabel(' Distance from center of cluster (cm)');
xlabel(' Flight Phase');


figure(); hold on;
plot(Loca(g2),zscore(smooth(meanC(g2),30,'lowess')),'LineWidth',3);
plot(Loca(g2),zscore(smooth(meanF(g2),30,'lowess')),'LineWidth',3);
legend('Calcium consistancy','Flight Variability');
ylabel('normalized units (z)');
xlabel(' Flight Phase');
title('Code Stability vs Flight Variability ');

