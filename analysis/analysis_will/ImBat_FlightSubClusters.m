function [output] = ImBat_FlightSubClusters
% check if there is significant subclustering of the flights if you first
% clutser the calcium timeseries, then check if there is significant
% difference in the flight clusters

% run in processed folder


% initiate vars

pvalTotal = [];


DIR = pwd;

files = dir(pwd);
files(ismember( {files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.


for i = 1 : length(subFolders)
    load([subFolders(i).name, '/CellReg_files/ROI_Data/Saved_Data/Aligned_Data.mat'])
    load([subFolders(i).name, '/CellReg_files/ROI_Data/cellRegistered.mat'])
    
    
    % Load Data
    % load (flightPaths, cell_registered_struct);
    
    % I. Bat Specific Details
    mnFlights_per_day(i,:)  = size(flightPaths.id,1)/max(flightPaths.day);
    mn_sessions(i,:) = max(flightPaths.day);
    
    
    % II. Across-Bat Summary Details
    
    % 1. Calculate tracked ROIs
    G = cell_registered_struct.cell_to_index_map;
    G(G>0) = 1;
    
    
    G2 = sum(G');
    idx2use = find(G2>1);
    
    for ii = 1:3
        for ix = 1:size(idx2use,2);
            [pval(ii,ix)]= ImBat_ClusterCalciumVar(FlightAlignedROI{ii},idx2use(ix));
            close all;
        end
    end
    
g = (pval<0.05/3); % index of significant one
g2 = sum(g);
g2(g2>0) = 1; % turn into sigificanct binary
    
    output.pval_per_bat{i} = pval;
    output.percent_Signif{i} = sum(g2)/size(g2,2)*100;
    
    pvalTotal = cat(2,pvalTotal,g2);
    clear pval;
end

output.pvalTotal = pvalTotal;
% report the % significantly sub-clustered
sum(output.pvalTotal)/size(output.pvalTotal,2)*100
