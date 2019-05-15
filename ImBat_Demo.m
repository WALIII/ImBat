function [out] = ImBat_Demo

% Run in extracted folder w/ tifs


% Align Video and save flight trajectories

disp('Aligning Flight trajectories');
[out] = ImBat_alignTimeStamps(audio,video,TS,Markers);

% load in audio:


% Load in all vid
