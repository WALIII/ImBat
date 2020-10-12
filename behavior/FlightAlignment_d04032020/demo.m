
% Behavioral ALignment Demo (run in FlightAlignment folder) 
% WAL3 04/03/2020

% Load example data 
out1 = load('example_data/alignment/Gal_200311_fly-1_Alignment.mat','out');
date1 = '200311';
Location_01 = out1.out.Location;

out2 = load('example_data/alignment/Gal_200315_fly-1_Alignment.mat','out');
Location_02 = out2.out.Location;
date2 = '200315';

% Load master Tracking File
MasterTrack = load('Master_Tracking_File.mat','out')


% Run alignment
Location_01_adjusted = ImBat_Align_Tracking(Location_01,MasterTrack.out,date1);
Location_02_adjusted = ImBat_Align_Tracking(Location_02,MasterTrack.out,date2);

% Plot unadjusted flights
figure();
hold on;
plot(Location_01(:,1),Location_01(:,2));
plot(Location_02(:,1),Location_02(:,2));
title('Unaligned data');


% Plot adjusted flights
figure();
hold on;
plot(Location_01_adjusted(:,1),Location_01_adjusted(:,2));
plot(Location_02_adjusted(:,1),Location_02_adjusted(:,2));
title('Aligned data');
