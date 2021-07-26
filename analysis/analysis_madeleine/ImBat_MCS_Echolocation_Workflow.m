% Workflow for extracting echolocations from microphone data and aligning
% it to flights.

clear all;

%% Variables to change:
thedate = '200518';
load(strcat('../extracted/Gen_',thedate,'_fly-1_track.mat'));
% CHANGE THIS TO MATCH THE FLIGHT DATA
%load(strcat('../extracted/Gen_',thedate,'_fly-1_extraction/processed_2020_12_18__1231/AV_data.mat'));
load(strcat('../extracted/processed_2021_01_28__0917/AV_data.mat'));

%% DO NOT CHANGE
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
    [echolocation_peaks,echolocation_vector_DS] = ImBat_MCS_extract_echolocations(audioConCat,out);

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

% Group the flights and in the flightpaths struct include the microphone
% data and echolocation data
clear ROI_Data
ROI_Data{1}.date = thedate;
ROI_Data{1}.Alignment = load('Alignment.mat');
ROI_Data{1, 1}.ROIs.results.metadata.cnmfe.Fs = 192000;
cluster_to_plot=1;
flightPaths34 = ImBat_MCS_GroupFlights(ROI_Data,cluster_to_plot,audioConCat,All_echolocation_vector_DS,'dist',2.1);

ImBat_MCS_
