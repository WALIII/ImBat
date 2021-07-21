% Workflow for extracting echolocations from flights

% Concatinate the microphone data and save it
% audio_concat on MCS_PC in lab
clear all; 
audiofile = 2;
thedate = '200518';

load(strcat('audioConCat_',num2str(audiofile),'.mat'));
load('ttlConCat.mat');
load(strcat('../extracted/Gen_',thedate,'_fly-1_track.mat'));
% CHANGE THIS
%load(strcat('../extracted/Gen_',thedate,'_fly-1_extraction/processed_2020_12_18__1231/AV_data.mat'));
load(strcat('../extracted/processed_2021_01_28__0917/AV_data.mat'));
% Re-create the alignment.mat file
% Manually load in
%   AV_data.mat
%   Alignment.mat
%   ...track.mat
%   audioConCat.mat
%   ttlConCat.mat

[echolocation_peaks] = ImBat_MCS_extract_echolocations(audioConCat);

TS = AnalogSignals;
[out, metrics] = ImBat_MCS_alignTimeStamps(audio,video,TS,Markers,audioConCat,ttlConCat,echolocation_peaks);

% Group the flights and in the flightpaths struct include the microphone
% data
clear ROI_Data
ROI_Data{1}.date = thedate;
ROI_Data{1}.Alignment = load('Alignment.mat');
ROI_Data{1, 1}.ROIs.results.metadata.cnmfe.Fs = 192000;
cluster_to_plot=1;
flightPaths34 = ImBat_MCS_GroupFlights(ROI_Data,cluster_to_plot,audioConCat,'dist',2.1);
temp_flight_echos = flightPaths34.echos;
temp_flight_echos_vector = flightPaths34.EcholocationIdx;
save(strcat('flight_echos_',num2str(audiofile)),'temp_flight_echos');
save(strcat('flight_echoVect_',num2str(audiofile)),'temp_flight_echos_vector');

ImBat_MCS_Echolocation_Analysis(flightPaths34)

num_clust = size(flightPaths34.clusterIndex,2);
% Create vector of all flight cluster IDs 
true_times = flightPaths34.AllFlightsMasterTime(flightPaths34.flight_starts_idx(:));
[tts,tts_idx]  = sort(true_times(:));
c_s_34 = flightPaths34.id(tts_idx);
