function [ROI_Data] = ImBat_MultiDayAnalysis
% First-pass Multi-day analysis function 

% WAL3
% d190609

% User input: processing steps to run ( t/f)
ST1 = 1; % extract data
ST2 = 1; % multiday alignment
ST3 = 1; % Across day analysis
    ST3_1 = 1; % Max projection overlay 
    ST3_2 = 1; % Place Cell Overlay
    ST3_3 = 1; %
    ST3_5 = 1;

% local Directory:
LD = '/Users/ARGO/Documents/DATA/Bat_Data';

% Cell array containing strings of the days you want to look at
days = {'190529','190530','190604'};

% Run in the 'processed' folder containing the days
DIR = pwd;

%% Go through all days and grab the data
if ST1 == 1;
for i = 1: size(days,2)
    % Save the date:
    ROI_Data{i}.date = days{i};
    
    % create index into correct subfolder
    folder = [days{i},'/','Zack_',days{i},'_fly-1_extraction/processed'];
    folder_rest = [days{i},'/','Zack_',days{i},'_rest1-1_extraction/processed'];
    
    disp(['entering into folder for day:  ', days{i}]);
    cd(folder) % index into folder 
    Alignment = load('Alignment.mat'); % get alignment data
    ROIs = load('Motion_corrected_Data_DS_results.mat'); % Get ROI data
   
    % Consolidate ROI and Flight data
    ROI_Data{i}.Alignment = Alignment;
    ROI_Data{i}.ROIs = ROIs;
    clear Alignment ROIs

    % Make Max projections ( Flights )
    disp( 'Getting Max Projection for Flights...');
    load('Motion_corrected_Data','Y');
    [MaxProj_flights, ~] = ImBat_Dff(Y);
    ROI_Data{i}.MaxProj_flights = MaxProj_flights;
    clear Y MaxProj_flights;
    
    % Make Max projections ( Rest )
    cd(DIR); cd(folder_rest);
    load('Motion_corrected_Data', 'Y');
    disp( 'Getting Max Projection for rest...');
    [MaxProj_rest, ~] = ImBat_Dff(Y);
    ROI_Data{i}.MaxProj_rest = MaxProj_rest;
    clear Y MaxProj_rest;

        cd(DIR)
end

% Save data locally somewhere...
cd(LD);
for i = 1:size(ROI_Data,2);
    if i ==1;
        filename = ['ROI_Data_',days{1}]
    else
        filename = [filename,'_',days{i}];
    end
end
    
save([filename,'.mat'],'ROI_Data')

end

%% Cell Sort to get 'similar ROIs'
if ST2 ==1;
    
end


%% Plot place cell-ness of the top most consistant cells
if ST3 == 1
 % Basic alignment:
 A = ROI_Data{2}.MaxProj_rest;
 B = ROI_Data{2}.MaxProj_flights;
[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = 0.009;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;
A = imregister(A, B, 'rigid', optimizer, metric);

[RGB1 RGB2] = CaBMI_XMASS(A,B,A);
figure(); imagesc(squeeze(RGB1))
 
end







