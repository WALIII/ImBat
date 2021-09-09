

% create random flights
%[FlightAlignedROI_rando] = ImBat_Align_FC(CombinedROI,flightPaths,1);

ds_factor = 10;

% get aligned data
[S2save,F2save,other] =  ImBat_Spikes(FlightAlignedROI);

% train on distance:
Y = other.T_dist;
Y = downsample(Y,ds_factor);

X = downsample(S2save,ds_factor);
X = round(X);

%X2 = X(end-3000:end,:);
%Y2 = Y(end-3000:end);

%X(end-3000:end,:) = [];
%Y(end-3000:end) = [];


Mdl = fitcnb(X,Y,'Distribution','kernal');
isGenRate = resubLoss(Mdl,'LossFun','ClassifErr');

label = predict(Mdl,X);;
figure();
hold on;
plot((label))
plot((Y-10))


% label2 = predict(Mdl,X2);
% figure();
% hold on;
% plot((label2))
% plot((Y2-10))
% title('witheld data');
% 


% %%%%
% 
% % apply now to random flights
% [FlightAlignedROI_rando] = ImBat_Align_FC(CombinedROI,flightPaths,1);
% FlightAlignedROI_rando1{1} = FlightAlignedROI_rando;
% [S2save_test,F2save_test,other_test] =  ImBat_Spikes(FlightAlignedROI_rando1);
% 
% 
% % predict position in flight in witheld data
% 
% X2 = downsample(S2save_test,ds_factor);
% X2 = round(X2);
% 
% % test on distance:
% Y2 = other_test.T_dist;
% Y2 = downsample(Y2,ds_factor);
% 
% 
% label2 = predict(Mdl,X2);
% figure();
% hold on;
% plot((label2))
% plot((Y2-10))
% title('witheld data');
% 
