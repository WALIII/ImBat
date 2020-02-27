function extract_tracking_data(varargin)

homeDir = pwd;
if nargin==2
    entry = 1;
else
    entry = 0;
end
% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'homedir'
            homeDir=varargin{i+1};
    end
end

if entry ==1
    c3dList = dir([homeDir filesep 'raw_tracking' filesep  '*.c3d']);
else
    c3dList = dir([homeDir filesep '*.c3d']);
end

%run through list of c3d files in that directory, convert to mat, and save
%to processed directory
for i = 1:length(c3dList)
    fileName = extractBefore(c3dList(i).name,'.c3d');  %'-Bat_Cluster.c3d'); %Bat
    dateSesh = datestr(datetime(c3dList(i).date),'yymmdd');
    batName = extractBefore(fileName,['_' dateSesh]);
    %sessionNum = fileName(end);
    %make processed directory
    if strcmp(batName,'setup')
        processed_dir = [homeDir filesep 'extraction_processing' filesep 'extracted_trackingCalibration' filesep];
    else
        processed_dir = [homeDir filesep 'extraction_processing' filesep 'extracted_scopeTracking' filesep batName(1:2) dateSesh filesep];
    end
    if ~isdir(processed_dir)
        mkdir(processed_dir);
    end
    if entry == 1
    [Markers,VideoFrameRate,AnalogSignals,AnalogFrameRate,Event,ParameterGroup,CameraInfo,ResidualError]=readC3D_analog([homeDir filesep 'raw_tracking' filesep c3dList(i).name]); %convert file
    else
    [Markers,VideoFrameRate,AnalogSignals,AnalogFrameRate,Event,ParameterGroup,CameraInfo,ResidualError]=readC3D_analog([homeDir filesep c3dList(i).name]); %convert file
    end
    %plot ttl impulses to check they are linear and not missing ttl
    event_ttls = AnalogSignals(:,2);
    [R,LT,UT,LL,UL] = risetime(event_ttls,VideoFrameRate);
    figure
    plot(LT)
    title(fileName)
    %save new mat file in both original directory and copied directory for processing
    save([processed_dir fileName '_track' '.mat'],'AnalogFrameRate','AnalogSignals','Markers','VideoFrameRate');
end