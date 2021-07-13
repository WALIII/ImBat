function out2 = ImBat_ROI_Behav_Correlation(FlightAlignedROI_combined)

% find if there is a correlation to the ROI and behavioral variance
close all
counter = 1;
min_flight_number = 30; % min number of flights needed
min_prom = 0.25;
disp(min_prom);

% NOTE: there may be more than one peak per ROI!

% to get data, run:
for ii = 1:size(FlightAlignedROI_combined{1}.C,1) % first 10 cells:
    out =  ImBat_analysis_11112020(FlightAlignedROI_combined,ii);
    
    % plot mean of ROI data
    % Zscore data:
    AllColZ = (out.data2.AllCol);
    %remove zeros
    xtemp = mean(AllColZ);
    clear h1 h2 Bpks Blocs h1 h2
    h1 = find(xtemp ==0);
    h2 = find(xtemp ~=0);
    AllColZ(:,h1) = [];
    clear x
    if size(AllColZ,2)<min_flight_number % remove if not enough flights
        x = zeros(size(AllColZ,1),size(AllColZ,2));
        x = mean(x,2);
    else
         x = mean(AllColZ,2);
    end
    
    
    % get peaks
    [Bpks,Blocs] = findpeaks(x,'MinPeakProminence',min_prom,'MinPeakDistance',120);
     clear x AllColZ
    % only use peaks in the flight pad ( exclude reward...
    % Blocs(Blocs<250) = [];
    % Blocs(Blocs>1200) = [];
    flightInterval = mean(FlightAlignedROI_combined{1}.FlightLength);
    stTime = 550;
    stpTime = stTime+flightInterval;
    Blocs(Blocs<stTime) = [];
    Blocs(Blocs>stpTime) = [];
    % Blocs(Blocs>1050) = [];
    numPeaks(ii) = size(Blocs,1);
    
    % figure();
    % hold on;
    % plot(x);
    %
    % plot(Blocs,ones(size(Blocs))+1,'*r');
    
    % now, get the distributions for 1. flights 2. calcium;
    if isempty(Blocs)
    else
        for iii = 1: size(Blocs,1);
            calPeaks = (mean(out.data2.AllCol(Blocs(iii)-20:Blocs(iii)+20,:)));
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
            
            % try to center on the field center...
            idx2use = find(calPeaks>prctile(calPeaks,90));
            %idx2use = h2;

            % get  position to mean point
            meanPoint(1) = median(flPeaksX(idx2use));
            meanPoint(2) = median(flPeaksY(idx2use));
            meanPoint(3) = median(flPeaksZ(idx2use));
            
            for i = 1: size(flPeaksZ,1)
                currPoint(1) = flPeaksX(i);
                currPoint(2) = flPeaksY(i);
                currPoint(3) = flPeaksZ(i);
                D(i) = vecnorm(currPoint - meanPoint, 2, 2);
            end
            
            % figure(); histogram(D);
            
            
            out2.D{counter} = D; % flight variance
            out2.calPeaks{counter} = calPeaks; % flight Peak
            out2.calPeakLoc{counter} = calPeakLoc; % flight Peak Location
            out2.calROI_ID{counter} = ii; % ROI_ID- possibly more than one peak per ROI
            counter = counter+1;
        end
    end
    
end
out2.numPeaks = numPeaks;

%% Plotting: ( will also work outside of function)
try
figure();
counter = 1;
for i = 1: size(out2.D,2)
    
    % remove zeros
    xx = out2.calPeaks{i};
    h = find(xx ==0);
    F2use = out2.D{i};
    C2use = out2.calPeaks{i};
    F2use(h) = [];
    C2use(h) = [];
    
    if size(F2use,2)>1;
    % normalize calcium from 0-1
     %C2use = [min_prom C2use]; % concat a zero to scale from 0->max
    C2use_raw = C2use;
    C2use = (C2use);
    
        %C2use(:,1) = []; % remove that zero..
    % ..now take the mean, and variance...
    
    varF(counter) = var(F2use);
    varC(counter) = std(C2use);
    varCr(counter) = std(C2use_raw);
    
    meanF(counter) = mean(F2use);%mean(zscore(F2use)-min(zscore(F2use)));
    meanC(counter) = mean(C2use);
    meanCr(counter) = mean(C2use_raw);
    FanoFactor(counter) = varF(counter)/meanF(counter);

    Loca(counter) = out2.calPeakLoc{i};
    %%% sig value:
    X1 = F2use';
    y = (C2use');
    % y(X1>300) = [];
    % X1(X1>300) = [];
    x1 = ones(size(X1,1),1);
    X = [x1 X1];    % Includes column of ones
    
    [~,~,~,~,stats] = regress(y,X);
    stats2keep(counter) = stats(3);


    counter = counter+1;
end
end
Loca = mat2gray(Loca);
%meanC = mat2gray(meanC);
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
stats
X = [X1];
mdl = fitlm(X,y);
figure();
plotAdded(mdl)
title('Calcium consistancy vs Flight Variance');
ylabel(' Mean, Normalized Ca+ activity');
xlabel('mean Flight distance (cm);');



% More plotting..... :

figure();
plot(meanF,(meanC),'o');
title('Mean Calcium intensity vs aveage flight distance from the mean flight');
ylabel(' Mean, Normalized Ca++ activity');
xlabel(' average euclidian distance from the mean flight in that cluster');



figure();
subplot(2,1,1);

hold on;
plot(Loca,(meanC),'o');
[g1 g2] = sort(Loca);
plot(Loca(g2),(smooth(meanC(g2),30,'lowess')),'LineWidth',3);

title('Mean Calcium intensity vs Flight Phase');
ylabel(' Mean, Normalized Ca++ activity');
xlabel(' Flight Phase');




subplot(2,1,2);
hold on;
plot(Loca,(meanF),'o');
plot(Loca(g2),(smooth(meanF(g2),30,'lowess')),'LineWidth',3);
title('Mean Flight distance from center of cluster vs Flight Phase');
ylabel(' Distance from center of cluster (cm)');
xlabel(' Flight Phase');


figure(); hold on;
plot(Loca(g2),zscore(smooth(meanC(g2),50,'rlowess')),'LineWidth',3);
plot(Loca(g2),zscore(smooth(meanF(g2),50,'rlowess')),'LineWidth',3);
legend('Calcium consistancy','Flight Variability');
ylabel('normalized units (z)');
xlabel(' Flight Phase');
title('Code Stability vs Flight Variability ');
catch
   out2 = [] 
end

try
out2.meanF = meanF;
out2.meanC = meanC;
out2.stats2keep = stats2keep;
catch
    out2.meanF = [];
out2.meanC = [];
out2.stats2keep = [];
end
