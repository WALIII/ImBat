function [R, R3, M] = ImBat_Flight_Model_AcrossDays(flightPaths1,day2use);


% model to test remapping in light/dark

% assumptions: HPC cares about a mix of things
% Location
% Head direction
% boundries
% time ( from takeoff)

% Default Model Paramaters:
sizeOfField = 500; % in mm
fireProb = 5; % 1 is highest 120 = 1 per s
numItterations = 2000;
numSpontanSpikes = 100; % max number of spontanious spikes
degree2use = 180;



% Plottinf things:
close all
alpha2use = 0.3;
col = lines(numItterations);
counter = 1;
Smarker = 100; % marker size


% get first and second third of data:
FF3 = flightPaths1.flight_starts_idx;
dc1 = flightPaths1.day;

% % randomize index
% fakeInd = randperm((size(dc1,1)));
% dc1 = dc1(fakeInd);


% for all flights
out.Light_all.FirstLight = find(dc1==1);
% out.Light_all.Darkness = find(dc>20 & dc<40);
out.Light_all.LastLight = find(dc1==day2use);

% temp = find(dc>40);
% temp2 = find(dc>4 & dc<20);
% temp = [temp' temp2'];
% temp = temp';
% out.Light_all.FirstLight = temp(1:2:end);
% out.Light_all.LastLight = temp(2:2:end);
% out.Light_all.Darkness = find(dc>20 & dc<40);



%figure(1);
hold on;
A = flightPaths1.tracjectoriesRaw*1000;
Flight1_concat = [];
Flight2_concat = [];

% Pre-Light
%figure(1);
subplot(1,2,1); hold on;
Ind2use = out.Light_all.FirstLight;
for iii = 1:length(Ind2use)
    bound = flightPaths1.flight_starts_idx(Ind2use(iii)):flightPaths1.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights
    Flight1_concat = cat(2,A(:,bound),Flight1_concat);
end
axis off
title('Lights ON');


% Light again
subplot(1,2,2); hold on;
Ind2use = out.Light_all.LastLight;
for iii = 1:length(Ind2use)
    bound = flightPaths1.flight_starts_idx(Ind2use(iii)):flightPaths1.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights
    Flight2_concat = cat(2,A(:,bound),Flight2_concat);
end
title('Lights On Again');


axis off

% now, generate a few 'place cells'

% %col = cat(2,col,ones(length(col),1)*0.1);
% hold on;
% subplot(1,2,1)
% plot3(Flight_Light_concat(1,:),Flight_Light_concat(2,:),Flight_Light_concat(3,:),'color',[0 0 0 0.1])
% hold on;
% subplot(1,2,2)
% plot3(Flight_Dark_concat(1,:),Flight_Dark_concat(2,:),Flight_Dark_concat(3,:),'color',[0 0 0 0.1])
% hold on;
for i = 1:numItterations
    clear AngleRange SpontaniousRate
    for xyz = 1:3;
        loc_temp= randi(5000,1)-2500;
        PlaceCells(i,:,xyz) = loc_temp;
    end
    temp = (randi(360)-180);
    AngleRange = temp:(temp+degree2use);
    [aa bb] = find(AngleRange>360);
    AngleRange(aa) = AngleRange(aa)-360;
    SpontaniousRate = randi(numSpontanSpikes,1);
    figure(1); subplot(1,2,1);
    
    clear temp aa bb
    b = [];
    a = [];
    
    [a,b] = find(Flight1_concat(1,:)>PlaceCells(i,:,1) & Flight1_concat(1,:)<PlaceCells(i,:,1)+sizeOfField ...
        & Flight1_concat(2,:)>PlaceCells(i,:,2) & Flight1_concat(2,:)<PlaceCells(i,:,2)+sizeOfField ...
        & Flight1_concat(3,:)>PlaceCells(i,:,3) & Flight1_concat(3,:)<PlaceCells(i,:,3)+sizeOfField);
    
    % Add random points
    b = [b randi(length(Flight1_concat(1,:)),SpontaniousRate,1)'];
    b(b<21) = [];
    b(b>length(Flight1_concat)-21) = [];
    
    % now check heading angle:
    for ii = 1:length(b)
        u(1) = Flight1_concat(1,b(ii)-20);
        u(2) = Flight1_concat(1,b(ii)+20);
        v(1) = Flight1_concat(2,b(ii)-20);
        v(2) = Flight1_concat(2,b(ii)+20);
        
        winddir = atan2d(u,v);
        Angle(ii) = winddir(1);
    end
    
    if length(b)>0
        dots2plot = find(Angle<max(AngleRange) & Angle>min(AngleRange));
        dots2plot= 1:fireProb:length(dots2plot);
        
         scatter3(Flight1_concat(1,b(dots2plot)),Flight1_concat(2,b(dots2plot)),Flight1_concat(3,b(dots2plot)),ones(length(b(dots2plot)),1)*Smarker,'MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerFaceAlpha',alpha2use,'MarkerEdgeAlpha',alpha2use)
        PM1 = ImBat_2dHeatMap(Flight1_concat,b(dots2plot));
        
        clear dots2plot Angle
        
    end
    
    
    
    % Light 2
    b3 = [];
    a = [];
  %  figure(1); subplot(1,2,2);
    [a,b3] = find(Flight2_concat(1,:)>PlaceCells(i,:,1) & Flight2_concat(1,:)<PlaceCells(i,:,1)+sizeOfField ...
        & Flight2_concat(2,:)>PlaceCells(i,:,2) & Flight2_concat(2,:)<PlaceCells(i,:,2)+sizeOfField ...
        & Flight2_concat(3,:)>PlaceCells(i,:,3) & Flight2_concat(3,:)<PlaceCells(i,:,3)+sizeOfField);
    
    b3 = [b3 randi(length(Flight2_concat(1,:)),SpontaniousRate,1)'];
    b3(b3<21) = [];
    b3(b3>length(Flight2_concat)-21) = [];
    
    % now check heading angle:
    for ii = 1:length(b3)
        u(1) = Flight2_concat(1,b3(ii)-20);
        u(2) = Flight2_concat(1,b3(ii)+20);
        
        v(1) = Flight2_concat(2,b3(ii)-20);
        v(2) = Flight2_concat(2,b3(ii)+20);
        
        winddir = atan2d(u,v);
        Angle(ii) = winddir(1);
    end
    
    if length(b3)>0
        dots2plot = find(Angle<max(AngleRange) & Angle>min(AngleRange));
        dots2plot= 1:fireProb:length(dots2plot);
        
         scatter3(Flight2_concat(1,b3(dots2plot)),Flight2_concat(2,b3(dots2plot)),Flight2_concat(3,b3(dots2plot)),ones(length(b3(dots2plot)),1)*Smarker,'MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerFaceAlpha',alpha2use,'MarkerEdgeAlpha',alpha2use)
        
        PM3 = ImBat_2dHeatMap(Flight2_concat,b3(dots2plot));

        clear dots2plot Angle
        
    end
    
    
    %figure(2); 
    subplot(1,1,1); 
    imagesc(PM1);  
    subplot(1,3,3);
    imagesc(PM3); 
    

    %% Calculation
    
    
    if isempty(b)==0 || isempty(b3)==0 
        if (size(b,2)+size(b3,2))/2 > numSpontanSpikes*1.5% rate is less than spon rate...
        try
            
            R(counter) = corr2(PM1,PM3); % compare light conditions
            if counter>1
                R3(counter) = corr2(PM1,PM3_last);
            end
            PM3_last = PM3;

            disp(['mean(Light/Light): ',num2str(mean(R))]);
            M(1,counter) = max(PM1(:));
            M(2,counter) = max(PM3(:));
            counter = counter+1;
            
        catch
            disp('oo')
        end
        
        end
    end
    
    
    
    
    
    
    

end
% end

% % remove ones
 R(R<.01) = [];
% R3(R3<.01) = [];




idx = find(M(1,:)>median(M(1,:))*1 | M(2,:)>median(M(2,:))*1 );


figure();
hold on;
histogram(1-R(:),'Normalization','probability','FaceColor','c','BinWidth',0.05);
% histogram(1-R2(:),'Normalization','probability','FaceColor','k','BinWidth',0.05)
histogram(1-R3(:),'Normalization','probability','FaceColor','r','BinWidth',0.05)

title('Apparent remapping due to changes in behavior');
legend('Light v Light', 'Light v Dark','ID Shuffle')
xlabel('Remapping Index [ 1- r1]');
ylabel('Probability');

% plot CDF:
figure();  hold on; plot(sort(R));  plot(sort(R3),'r')
title('CDF of remapping due to changes in behavior');
legend('Light v Light', 'Light v Dark','ID Shuffle')
end


