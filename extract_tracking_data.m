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
 fileName = c3dList(i).name;
 batName = extractBefore(fileName,'_');
 sessionNum = extractBefore(fileName,'-');
 dateSesh = datestr(datetime(c3dList(i).date),'yymmdd');
[Markers,VideoFrameRate,AnalogSignals,AnalogFrameRate,Event,ParameterGroup,CameraInfo,ResidualError]=readC3D_analog([wd filesep fileName]); %convert file
%plot ttl impulses to check they are linear and not missing ttl
event_ttls = AnalogSignals(:,2);
[R,LT,UT,LL,UL] = risetime(event_ttls,VideoFrameRate);
figure
plot(LT)
title(fileName)
%save new mat file
save([processed_dir batName '_' dateSesh '_track_' sessionsNum '.mat'],'AnalogFrameRate','AnalogSignals','Markers','VideoFrameRate')
end