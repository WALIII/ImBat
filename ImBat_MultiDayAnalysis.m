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
saveFlag = 1;

% local Directory:
LD = 'Z:\users\tobias\flight\data_processed\topQualityData\analysis_done\za\';
%LD =  'D:\Bat_Data2\Zack\Processed\'
%LD = '/Users/ARGO/Documents/DATA/Bat_Data_ZuZu';
if saveFlag == 1
    %saveDir1 = '/Volumes/Tobias_flig/topQualityData/analysis_done/plots/';
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\ForTobias\plots\';
    %saveDir1 = '/Users/periscope/Desktop/analysis/flight/plots/';
    %saveDir1 = 'C:\Users\tobias\Desktop\analysis\plots\';
    if ~exist([saveDir1 datestr(now,'yymmdd') ])
        mkdir([saveDir1 datestr(now,'yymmdd')]);
    else
        disp('You have been working today...');
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep];
end
if ~exist(LD)
mkdir(LD);
end

% Cell array containing strings of the days you want to look at ( MANUAL)
% days = {'190528','190529','190530'};
% Get all folders in directory
files = dir(pwd);
files(ismember( {files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.
for k = 1 : length(subFolders)
	fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
    days{k} = subFolders(k).name;
end





% Run in the 'processed' folder containing the days
DIR = pwd;

% check into DIRs
cd(LD);
cd(DIR);

%% Go through all days and grab the data
if ST1 == 1;
for i = 1: size(days,2)
    % Save the date:
    ROI_Data{i}.date = days{i};
    
    % create index into correct subfolder ( get biggest flight and rest
   % Flight Flder
   
   cd([days{i} filesep 'extracted']);
   out = ImBat_GetBiggestFolder('fly');
    cd(out);
   processedDir = dir('processed_*');
    folder = [processedDir(end).name];
    ROI_Data{i}.folder = folder; % save this for indexing later...
    
    disp(['entering into folder for day:  ', days{i}]);
    cd(folder) % index into folder 
    Alignment = load('Alignment.mat'); % get alignment data
    ROIs = load('results.mat'); % Get ROI data
   
    % Consolidate ROI and Flight data
    ROI_Data{i}.Alignment = Alignment;
    ROI_Data{i}.ROIs = ROIs;
    clear Alignment ROIs

    % Make Max projections ( Flights )
    disp( 'Getting Max Projection for Flights...');
    load('Motion_corrected_Data_DS','Y','Ysiz');
    [MaxProj_flights, ~] = ImBat_Dff(Y);
    ROI_Data{i}.MaxProj_flights = MaxProj_flights;
    ROI_Data{i}.Ysiz_DS = Ysiz;
    clear Y Ysiz MaxProj_flights;
    
    
    % Make Max projections ( Rest )

        
try
 %cd(days{i});
 cd ..;
 cd ..;
 out = ImBat_GetBiggestFolder('rest1');  
cd(out);
   processedDir = dir('processed_*');
    folder_rest = [processedDir(end).name];
    
    %cd(DIR);
    cd(folder_rest);
    load('Motion_corrected_Data_DS', 'Y','Ysiz');
    disp( 'Getting Max Projection for rest...');
    [MaxProj_rest, ~] = ImBat_Dff(Y);
    ROI_Data{i}.MaxProj_rest = MaxProj_rest;
    ROI_Data{i}.Ysiz_DS = Ysiz;
    clear Y Ysiz MaxProj_rest;
catch
    disp('No rest data');
end

        cd(DIR)
end

% Save data locally somewhere...
cd(saveDir);
for i = 1:size(ROI_Data,2);
    if i ==1;
        filename = ['ROI_Data_',days{1}]
    else
        filename = [filename,'_',days{i}];
    end
end

save([filename, '.mat'],'ROI_Data')
%save([filename,'.mat'],'ROI_Data')
cd(LD)
end

%% Cell Sort to get 'similar ROIs'
if ST2 ==1;
    
end


%% Plot place cell-ness of the top most consistant cells
if ST3 == 1
%  % Basic alignment within day:
%  figure();
counter = 1;
for i = 1:size(days,2)-2;

%  clear A B C
 % Across day alignment ( flights)
 figure(); 

 A = ROI_Data{counter}.MaxProj_flights;
 B = ROI_Data{counter+1}.MaxProj_flights;
 C = ROI_Data{counter+2}.MaxProj_flights;
 
% A = imregister(A, B, 'rigid', optimizer, metric);
% C = imregister(C, B, 'rigid', optimizer, metric);
[RGB1 RGB2] = CaBMI_XMASS(A,B,C);
imagesc(squeeze(RGB1))

counter = counter+1;
end


% Condition Registration by saving data in the proper format
mkdir('CellReg_files');
cd('CellReg_files')

for i = 1:size(days,2)
    
    % get size of ROI field
    A2 = full(ROI_Data{i}.ROIs.results.A);
    for ii = 1:size(A2,2);
        A3(ii,:,:) = reshape(A2(:,ii),ROI_Data{i}.Ysiz_DS(1),ROI_Data{i}.Ysiz_DS(2)); 
    end
         A2 = A3(1:ii,:,:); % ROI masks in CellReg format...

    % save data:
    save([ROI_Data{i}.date,'.mat'],'A2');
    clear A2 A3
    
end


mkdir('ROI_Data');
save('ROI_Data/ROI_Data','ROI_Data');
% Register Cells by ROI masks:
CellReg

% Save indexes for similar ROIs, and load into ROI_Data ( for across-day
% indexing)



end







