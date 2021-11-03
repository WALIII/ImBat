% Workflow for extracting echolocations from microphone data and aligning
% it to flights.

%% Variables to change:
clear all;
%thedate = '200518';
%thedate = '200305';
thedate = '200518';

load(strcat('../extracted/Gen_',thedate,'_fly-1_track.mat'));
% CHANGE THIS TO MATCH THE FLIGHT DATA
%load(strcat('../extracted/Gen_',thedate,'_fly-1_extraction/processed_2020_12_18__1231/AV_data.mat'));
%Ge200518 - 
load(strcat('../extracted/processed_2021_01_28__0917/AV_data.mat'));
%Ge200519 - 
%load(strcat('../extracted/processed_2021_01_28__1843/AV_data.mat'));
%% Align and extract echolocations
All_echolocation_vector_DS = [];
for xx=1:5
    
    audiofile = xx;

    % Load in the audioConCat file and ttlConCat file
    load(strcat('audioConCat_',num2str(audiofile),'.mat'));
    load('ttlConCat.mat');

    % Align the Cortex TTL (TS), Framegrabber TTL (audio), Microphone TTL (ttl)
    TS = AnalogSignals;
    [out, metrics] = ImBat_MCS_alignTimeStamps(audio,video,TS,Markers,audioConCat,ttlConCat);

    % Extract echolocation peaks from the audioConCat data
    tic;
    [echolocation_peaks,echolocation_vector_DS] = ImBat_MCS_extract_echolocations(audioConCat,out);
    toc
    
    disp("Adding");
    % Add this microphone's echolocation indexes to a joint vector
    All_echolocation_vector_DS = [All_echolocation_vector_DS,echolocation_vector_DS];

    % Plot things to make sure the mic/ttl alignment looks correct
    figure(); hold on;
    title(strcat("Microphone #",num2str(xx)));
    plot(metrics.E_TS_tv_offset,metrics.TS_infer+0.02);
    plot(metrics.E_ttl_tv_offset,metrics.ttl_infer);
    plot(out.Microphone_Time,downsample(audioConCat,1600));
    plot(out.Microphone_Time,echolocation_vector_DS,'*r');
    clear audioConCat out metrics;
end

%% Group and Plot flights and echolocations
% Group the flights and in the flightpaths struct include the microphone
% data and echolocation data
clear ROI_Data
load('audioConCat_2.mat');
ROI_Data{1}.date = thedate;
ROI_Data{1}.Alignment = load('Alignment.mat');
try
    load('All_echolocation_vector_DS.mat');
catch
    disp("Run Code Above, no All_echolocations Vector");
end
ROI_Data{1, 1}.ROIs.results.metadata.cnmfe.Fs = 192000;
cluster_to_plot=1;
flightPaths34 = ImBat_MCS_GroupFlights(ROI_Data,cluster_to_plot,audioConCat,All_echolocation_vector_DS,'dist',2);

%% Questions about echolocation occurrance 

% Get the difference in echolocations across light and "dark" conditions
[dark_flights_idx,lite_flights_idx] = ImBat_MCS_light_conditions(flightPaths34);
% Plot the echolocations as a raster for all flights
ImBat_MCS_RMS_light_conditions(flightPaths34,dark_flights_idx,lite_flights_idx);
% Calculate the mean amplitde of the peaks summed across all flights
% normalized by the number of flights. 
cluster=8;
ImBat_MCS_echolocation_stereotypy(flightPaths34,cluster);
% For light-only sessions 
ImBat_MCS_conditions(flightPaths34,1);
ImBat_MCS_RMS(flightPaths34);
