

% create random flights
[FlightAlignedROI_rando] = ImBat_Align_FC(CombinedROI,flightPaths,1);

% random flights
for ii = 1:100
Flights = [];
Spikes =  [];

hold on;
bound2use = 1:1400;

cell2use = ii;

% get flight data, 
ClustFlight = FlightAlignedROI_rando.ClustFlight_withPads;

% upsample the calcium
CutCells = FlightAlignedROI_rando.S;
%CutCells = FlightAlignedROI{iii}.S;

for i = 1: size(CutCells,3);
    
    bound2use = 1:1400;
trial2use = i;
exampFlight = ClustFlight(bound2use,:,trial2use);
exampCell = CutCells(cell2use,:,trial2use);
exampCell = interp(exampCell,4);
exampCell = exampCell(1:size(exampFlight,1));

Flights = cat(1,Flights, exampFlight);
Spikes = cat(1,Spikes, exampCell');

end
% binarie spikes
Spikes =  smooth(zscore(Spikes),100);
% Spikes(Spikes>1) = 1;
% Spikes(Spikes<1) = 0;


S2save_train(:,ii) = Spikes;
F2save_train = Flights;
end

% structured flights
for ii = 1:100
Flights = [];
Spikes =  [];

hold on;
bound2use = 1:1400;

cell2use = ii;

% get flight data, 
ClustFlight = FlightAlignedROI{1}.ClustFlight_withPads;

% upsample the calcium
CutCells = FlightAlignedROI{1}.S;
%CutCells = FlightAlignedROI{iii}.S;

for i = 1: size(CutCells,3);
    
    bound2use = 1:1400;
trial2use = i;
exampFlight = ClustFlight(bound2use,:,trial2use);
exampCell = CutCells(cell2use,:,trial2use);
exampCell = interp(exampCell,4);
exampCell = exampCell(1:size(exampFlight,1));

Flights = cat(1,Flights, exampFlight);
Spikes = cat(1,Spikes, exampCell');

end
% binarie spikes
Spikes =  smooth(zscore(Spikes),100);
% Spikes(Spikes>1) = 1;
% Spikes(Spikes<1) = 0;


S2save_test(:,ii) = Spikes;
F2save_test = Flights;
end


% package data:
X = S2save_train;
for i = 1:size(X,2);
X(:,i) = smooth(X(:,i),10);
end
X = downsample(X,10);
X(X<0) = 0;
X= X*10;
X = round(X);


X2 = S2save_test;
for i = 1:size(X2,2);
X2(:,i) = smooth(X2(:,i),10);
end
X2 = downsample(X2,10);
X2(X2<0) = 0;
X2= X2*10;
X2 = round(X2);


% Decoder stuff:


% predict X and Y coordinates with 2 models 


%mdl_glm = fitglm(X,round(Y)+1,'Distribution','poisson');
% Get XY data:


xdata = F2save_train(1:end-1,1);
ydata = F2save_train(1:end-1,2);

xdata = downsample(smooth(xdata',10),10);
ydata = downsample(smooth(ydata',10),10);



% fit model
mdl_X = fitcnb(X,xdata,'Distribution','kernel');
mdl_Y = fitcnb(X,ydata,'Distribution','kernel');

% predict velocity
label_X = predict(mdl_X,X);
label_Y = predict(mdl_Y,X);

% withheld data
xdata_T = F2save_test(1:end-1,1);
ydata_T = F2save_test(1:end-1,2);

xdata_T = downsample(smooth(xdata_T',10),10);
ydata_T = downsample(smooth(ydata_T',10),10);

% predict based on witheld data
label_X2 = predict(mdl_X,X2);
label_Y2 = predict(mdl_Y,X2);


% Plot predictions
figure()
hold on;
plot(ydata,xdata,'b'); 

plot(label_Y,label_X,'r')




% % Model:
% % Y = smooth(F2save_train(:,1)',50);
% % Y = downsample(Y,50)/10;
% % Y = round(smooth(Y,10));
% Y = round(downsample(F2save_train(:,1),10)/10);
% X = S2save_train;
% for i = 1:size(X,2);
% X(:,i) = smooth(X(:,i),10);
% end
% X = downsample(X,10);
% X(X<0) = 0;
% X= X*10;
% X = round(X);
% 
% Mdl = fitcnb(X,Y,'Distribution','mn');
% isGenRate = resubLoss(Mdl,'LossFun','ClassifErr')

% isLabels1 = resubPredict(Mdl);
% ConfusionMat1 = confusionchart(Y,isLabels1);
% 
% % train Velocity
% x1 = F2save_train(2:end,1);
% x2 = F2save_train(1:end-1,1);
% y1 = F2save_train(2:end,2);
% y2 = F2save_train(1:end-1,2);
% 
% S=sqrt((x2-x1).^2+(y2-y1).^2);
% S(S>40) = 0;
% 
% Y = smooth(S',10);
% Y = downsample(Y,10);
% 
% 
% X = downsample(S2save_train,10);X(X<0) = 0;
% X= X*10;
% X = round(X);
% 
% Mdl = fitcnb(X,Y,'Distribution','poisson');
% isGenRate = resubLoss(Mdl,'LossFun','ClassifErr');
% 
% label = predict(Mdl,X);;
% figure();
% hold on;
% plot((label))
% plot((Y-10))
% 
% % test Velocity
% clear x1 y1 x2 y2 S;
% x1 = F2save_test(2:end,1);
% x2 = F2save_test(1:end-1,1);
% y1 = F2save_test(2:end,2);
% y2 = F2save_test(1:end-1,2);
% 
% S=sqrt((x2-x1).^2+(y2-y1).^2);
% S(S>40) = 0;
% 
% Y2 = smooth(S',10);
% Y2 = downsample(Y2,10);
% X2 = downsample(S2save_test,10);X2(X2<0) = 0;
% X2= X2*10;
% X2 = round(X2);
% 
% % predict velocity
% label2 = predict(Mdl,X2);
% figure();
% hold on;
% plot((label2))
% plot((Y2))
% 
% 
% label1 = predict(Mdl,X);
% figure();
% hold on;
% plot((label1))
% plot((Y))
% 
% 
% 
% Model_performance = corr(Y2,label2);
% 
% 
% 
% %%%%
% 
% 
% % 
% % 
% % % X(idx2rm,:) = [];
% % 
% % 
% % Y2 = smooth(F2save_test(:,1)',10);
% % Y2 = downsample(Y2,10);
% % Y2 = round(Y2/1);
% % idx2rm2 = find(abs(diff(Y2))==0);
% % Yy = Y2;
% % Yy2(idx2rm2) = [];
% % 
% % 
% % 
% % 
% % X2 = downsample(S2save_test,10);X2(X2<0) = 0;
% % X2= X2*10;
% % X2 = round(X2);
% % X2(idx2rm2,:) = [];
% % 
% % 
% % Y3 = cat(1,Y,Y2);
% % X3 = cat(1,X,X2);
% % 
% % label = predict(Mdl,X);
% % 
% % Mdl3 = fitcnb(X3,Y3,'Distribution','mn');
% % isGenRate = resubLoss(Mdl,'LossFun','ClassifErr')
% % 
% % % notPi = ~ismember(Y2, Y);
% % % ind2rm = find(notPi==1);
% % 
% % % Y2(ind2rm,:) = [];
% % % X2(ind2rm,:) = [];
% % oosGenRate = loss(Mdl,X2,(Y2));
% % 
% % label = predict(Mdl,X);
% % 
% % label3 = predict(Mdl3,X3);
% % 
% % 
% % figure();
% % hold on;
% % plot(smooth(label3,30))
% % plot(smooth(Y3,10))
% % 
% % 
% % 
% % figure();
% % hold on;
% % plot(smooth(label,30))
% % plot(smooth(Y,10))
% % 
% % figure();
% % hold on;
% % plot(smooth(label2,30))
% % plot(smooth(Y2,10))
% % 
% % %Y = round(Y/1);
% % % idx2rm = find(abs(diff(Y))==0);
% % % Yy = Y;
% % % Yy(idx2rm) = [];
% % 
% % % Y = Yy;
% % % clear Yy;
