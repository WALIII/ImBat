function extract_tracking_data(wd)
c3dList = dir([wd filesep '*.c3d']);
%make processed directory
processed_dir = [wd filesep 'processed' filesep];
if ~isdir(processed_dir)
    mkdir(processed_dir);
end


%run through list of c3d files in that directory, convert to mat, and save
%to processed directory
for i = 1:length(c3dList)
    fileName = extractBefore(c3dList(i).name,'-Bat_Cluster.c3d'); %Bat
    %fileName = extractBefore(c3dList(i).name,'-Bat_Cluster.c3d'); %Bat
    dateSesh = datestr(datetime(c3dList(i).date),'yymmdd');
    batName = extractBefore(fileName,['_' dateSesh]);
    %sessionNum = fileName(end);
    copy_dir = [extractBefore(wd,batName) 'processed' filesep batName filesep dateSesh filesep];
    if ~isdir(copy_dir)
        mkdir(copy_dir);
    end
    [Markers,VideoFrameRate,AnalogSignals,AnalogFrameRate,Event,ParameterGroup,CameraInfo,ResidualError]=readC3D_analog([wd filesep c3dList(i).name]); %convert file
    %plot ttl impulses to check they are linear and not missing ttl
    event_ttls = AnalogSignals(:,2);
    [R,LT,UT,LL,UL] = risetime(event_ttls,VideoFrameRate);
    figure
    plot(LT)
    title(fileName)
    %save new mat file in both original directory and copied directory for processing
    save([processed_dir fileName '_track' '.mat'],'AnalogFrameRate','AnalogSignals','Markers','VideoFrameRate');
    save([copy_dir fileName '_track' '.mat'],'AnalogFrameRate','AnalogSignals','Markers','VideoFrameRate');
end