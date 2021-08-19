

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

Mdl = fitcnb(X,Y,'Distribution','kernal');
isGenRate = resubLoss(Mdl,'LossFun','ClassifErr');

label = predict(Mdl,X);;
figure();
hold on;
plot((label))
plot((Y-10))


%%%%

% apply now to random flights
[FlightAlignedROI_rando] = ImBat_Align_FC(CombinedROI,flightPaths,1);
[S2save_test,F2save_test,other_test] =  ImBat_Spikes(FlightAlignedROI_rando);


% predict position in flight in witheld data

X2 = downsample(S2save_test,ds_factor);
X2 = round(X2);

% test on distance:
Y2 = other_test.T_dist;
Y2 = downsample(Y2,ds_factor);


label2 = predict(Mdl,X2);
figure();
hold on;
plot((label2))
plot((Y2-10))
title('witheld data');

