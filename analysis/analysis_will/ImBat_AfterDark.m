function [dark_cluster, out] = ImBat_AfterDark(flightPaths,varargin)
% Find flights that happened in the dark:

% Assumption is the lights are turned off exactly 20min, then back on at
% 40min
close all

% Get flightpaths of type:
FF = flightPaths.flight_starts_idx;
FlightPaths = 2;
pad2use = 100;

subplotting = 1; 

% User Inputs:
nparams=length(varargin);
if mod(nparams,2)>0
    error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
    switch lower(varargin{i})
        case 'flightpaths'
            FlightPaths=varargin{i+1};
    end
end
%figure(); plot(flightPaths.AllFlightsTime)

% plot all flights
figure();
histogram(flightPaths.AllFlightsTime(FF)/60,100);

% plot unclustered
for i = 1:size(FlightPaths,2)+1
    cluster = i;
    flightPaths.clusterIndex{cluster};
    FF2 = flightPaths.flight_starts_idx(flightPaths.clusterIndex{cluster});
    
    figure();
    hold on;
    histout{i} = histogram(flightPaths.AllFlightsTime(FF2)/60,'NumBins',60,'BinLimits',[1,60]);
    xlabel('time (min)');
    ylabel('count');
    title(['Cluster ',num2str(i)]);
    plot([20 20],[0 10],'--r')
    plot([40 40],[0 10],'--r')
    clear FF2
    xlim([0 60]);
end

% All flights:
figure();
FF2 = flightPaths.flight_starts_idx(:);
histout_all = histogram(flightPaths.AllFlightsTime(FF2)/60,'NumBins',60,'BinLimits',[1,60]);





% Export Data
for i = 1:size(FlightPaths,2)+1; % for cluster
    FF2 = flightPaths.flight_starts_idx(flightPaths.clusterIndex{i});
    dc = flightPaths.AllFlightsTime(FF2)/60;
    out.Light_cluster{i}.FirstLight = flightPaths.clusterIndex{i}(find(dc<20));
    out.Light_cluster{i}.Darkness = flightPaths.clusterIndex{i}(find(dc>20 & dc<40));
    out.Light_cluster{i}.LastLight = flightPaths.clusterIndex{i}(find(dc>40));
    
    dark_cluster{i} = find(dc>20 & dc<40); % would use out.Light_cluster{i}.Darkness' to index into the flight order... 
    clear FF2 dc;
end
FF3 = flightPaths.flight_starts_idx;
dc = flightPaths.AllFlightsTime(FF3)/60;

% for all flights
out.Light_all.FirstLight = find(dc<20);
out.Light_all.Darkness = find(dc>20 & dc<40);
out.Light_all.LastLight = find(dc>40);
out.HistData = histout;



% some basic plotting
%% Plot unclustered to top 3 Ratio
try
figure();
hold on;
%plot(medfilt1(histout{1}.Values./(histout{1}.Values+histout{2}.Values+histout{3}.Values),1,[],2))
plot(smooth(histout{1}.Values./(histout{1}.Values+histout{2}.Values+histout{3}.Values+histout{4}.Values),3),'--r')

xlabel('time (min)');
ylabel('proportion of unique flights');
title(['proportion of unique flights to top flights']);
plot([20 20],[0 1],'-b')
plot([40 40],[0 1],'-b')
catch
    disp('not enought flightpaths to plot stats');
end


%% Plot Flights

figure();
hold on;
A = flightPaths.tracjectoriesRaw*1000;

% Pre-Light
if subplotting ==1;
subplot(1,3,1); 
else
    figure();
end

hold on;
Ind2use = out.Light_all.FirstLight;
for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii))-pad2use:flightPaths.flight_ends_idx(Ind2use(iii))+pad2use;
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights
end
axis tight
axis off

title('Lights ON');
% Darkness
if subplotting ==1;
subplot(1,3,2); 
else
    figure();
end
hold on;
Ind2use = out.Light_all.Darkness;
for iii = 1:length(Ind2use)
    bound = flightPaths.flight_starts_idx(Ind2use(iii))-pad2use:flightPaths.flight_ends_idx(Ind2use(iii))+pad2use;
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights
end
title('Lights Off');
axis tight
axis off

% Last Light
if subplotting ==1;
subplot(1,3,3); 
else
    figure();
end
hold on;
Ind2use = out.Light_all.LastLight;
for iii = 1:length(Ind2use)
    bound = flightPaths.flight_starts_idx(Ind2use(iii))-pad2use:flightPaths.flight_ends_idx(Ind2use(iii))+pad2use;
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights
end
title('Lights Return on');
axis tight
axis off
%% Now, Plot the same figure, but overlay the clustered flights:

%% %% %% Plot Flights


figure();
col = {[1 0 0 0.2],[0 0 1 0.2],[0 1 0 0.2],[1 0 1 0.2],[1 1 0 0.2],[0 1 1 0.2]};

for i = 1:length(FlightPaths);
cluster = FlightPaths(i);
hold on;
A = flightPaths.tracjectoriesRaw*1000;

% Pre-Light
if subplotting ==1;
subplot(1,3,1); 
else
    figure(30);
end
hold on;
Ind2use = out.Light_all.FirstLight;
Ind2use2 = out.Light_cluster{cluster}.FirstLight  ;
axis off
if i ==1;
for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii))-pad2use:flightPaths.flight_ends_idx(Ind2use(iii))+pad2use;
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights

end
end
for ii = 1: length(Ind2use2);
    bound2 = flightPaths.flight_starts_idx(Ind2use2(ii))-pad2use:flightPaths.flight_ends_idx(Ind2use2(ii))+pad2use;
    plot2 =  plot3(A(1,bound2),A(2,bound2),A(3,bound2),'color',col{i},'LineWidth',2); % plot all flights
    Flights2save{i}{1}(:,:,ii) = A(:,bound2(1):bound2(1)+600);
end
title('Lights ON');
axis tight
axis off


% Darkness
if subplotting ==1;
subplot(1,3,2); 
else
    figure(31);
end
hold on;
Ind2use = out.Light_all.Darkness;
Ind2use2 = out.Light_cluster{cluster}.Darkness;
if i ==1;
for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii))-pad2use:flightPaths.flight_ends_idx(Ind2use(iii))+pad2use;
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights
end
end
for ii = 1: length(Ind2use2);
    bound2 = flightPaths.flight_starts_idx(Ind2use2(ii))-pad2use:flightPaths.flight_ends_idx(Ind2use2(ii))+pad2use;
    plot2 =  plot3(A(1,bound2),A(2,bound2),A(3,bound2),'color',col{i},'LineWidth',2); % plot all flights
    Flights2save{i}{2}(:,:,ii) = A(:,bound2(1):bound2(1)+600);
end
title('Lights Off');
axis tight
axis off


% Last Light
if subplotting ==1;
subplot(1,3,3); 
else
    figure(33);
end
hold on;
Ind2use = out.Light_all.LastLight;
Ind2use2 = out.Light_cluster{cluster}.LastLight;
if i ==1;
for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii))-pad2use:flightPaths.flight_ends_idx(Ind2use(iii))+pad2use;
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights
end
end
for ii = 1: length(Ind2use2);
    bound2 = flightPaths.flight_starts_idx(Ind2use2(ii))-pad2use:flightPaths.flight_ends_idx(Ind2use2(ii))+pad2use;
    plot2 =  plot3(A(1,bound2),A(2,bound2),A(3,bound2),'color',col{i},'LineWidth',2); % plot all flights
    Flights2save{i}{3}(:,:,ii) = A(:,bound2(1):bound2(1)+600);
end
title('Lights Return on');
end
axis tight
axis off


out.Flights2save = Flights2save; %data on flights

% reate custom color:
col2use(1,:) = [0 0 0];
col2use(2,:) = [1 0 0];
col2use(3,:) = [0 0 1];
col2use(4,:) = [1 0 1];
col2use(5,:) = [0 1 0];
col2use(6,:) = [0 1 1];
col2use(7,:) = [1 1 0];



clear a
for i = 1:size(FlightPaths,2)+1;
a(:,i) = histout{i}.Values/max(flightPaths.day);
end
figure();
hold on;
b = bar(a,'stacked');
    plot([20 20],[0 30/max(flightPaths.day)],'--k','LineWidth',2)
    plot([40 40],[0 30/max(flightPaths.day)],'--k','LineWidth',2)
title(' Number of flights vs time of day')
xlabel('time in session ( min)');
ylabel('avg flights/min');
legend('Unique','FlightPath 1','FlightPath 2','FlightPath 3','FlightPath 4');


for K = 1 : length(b);
    if K ==1;
       b(K).FaceAlpha = 0.5;
       b(K).EdgeAlpha = 0.5;
    else
       b(K).FaceAlpha = 0.8;
       b(K).EdgeAlpha = 0.8;
    b(K).FaceColor = col2use(K,:).'; 
    b(K).EdgeColor = col2use(K,:).'; 
    end
end


% Plot overlays:
col2(1,:) = [0 1 1];
col2(2,:) = [1 0 0];
col2(3,:) = [0 0 1];

figure();
hold on;
for ii = 1:3; 
    for iii = 1:3;
        subplot(1,3,iii);
        hold on;
adata =  squeeze(Flights2save{1}{ii}(iii,:,:))'+iii;
L = size(adata,2);
se = std(adata);%sqrt(length(adata));
mn = mean(adata);
mn = smooth(mn,1)';
h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col2(ii,:)); alpha(0.5);
plot(mn,'Color',col2(ii,:));
    end
end

%% Stability across phases:
% concat
% correlation
counter = 1;
for flight2use = 1:length(FlightPaths);
    clear Mf
Mf = cat(3,Flights2save{flight2use}{1}(:,:,:),Flights2save{flight2use}{2}(:,:,:),Flights2save{flight2use}{3}(:,:,:));
MF = squeeze(mean(Mf,3));

for LDL = 1:3; % 1: light, 2: dark, 3:light2
for i = 1: size(Flights2save{flight2use}{LDL}(:,:,:),3)
    try
MC_fX = corrcoef(squeeze(Flights2save{flight2use}{LDL}(1,:,i)),MF(1,:));
MC_fY = corrcoef(squeeze(Flights2save{flight2use}{LDL}(2,:,i)),MF(2,:));
MC_fZ = corrcoef(squeeze(Flights2save{flight2use}{LDL}(3,:,i)),MF(3,:));
MC_fX(2);
MC_f_all(counter,LDL) = (MC_fX(2)+MC_fY(2)+ MC_fZ(2))/3;
counter = counter+1;
    catch
        disp(' no flights in this cluster');
    end
end
end
end

figure();
MC_f_all(MC_f_all == 0) = NaN;
plotSpread(MC_f_all);


    out.PlotSpread = MC_f_all;
    

