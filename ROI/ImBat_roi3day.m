[roi3day] = ImBat_roi3day(videoData,flightPaths,alignment,varargin)

batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0;

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'batname'
            batName=varargin{i+1};
        case 'datesesh'
            dateSesh = varargin{i+1};
        case 'sessiontype'
            sessionType = varargin{i+1};
                    case 'loadflag'
            loadFlag = varargin{i+1};
    end
end

% ImBat_Dff
scaling = 4;

%labels for loading and saving data if running independent fromImBat_analyze
if loadFlag == 1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
    processedFolders = dir('processed*');
    processedNewest = sort({processedFolders(end).name});
    processedNewest = char(processedNewest);
    load([pwd processed_Newest '/Motion_corrected_Data_DS.mat']);
    analysisFolders = dir('analysis*');
    analysisNewest = sort({processedFolders(end).name});
    analysisNewest = char(processedNewest);
    load([pwd processed_Newest '/Motion_corrected_Data_DS.mat']);
    disp('')
end