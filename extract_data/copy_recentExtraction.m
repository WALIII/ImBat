dirBats = dir('z*');
%dirZbats = dir('z*');
homeDir = pwd;
saveDir = 'F:\flight\topQualityData\copydir';
 for bat_i = 1%:length(dirBats)
   cd(dirBats(bat_i).name);
   batDir = pwd;
   dirBatDays = dir([dirBats(bat_i).name(1:2) '*']);
   for day_i = 1:length(dirBatDays)
        cd([dirBatDays(day_i).name filesep 'extracted']);
        dayDir = pwd;
        
        %find and copy over tracking files
        dirTracking = dir('*_track.mat');
        for track_i = 1:length(dirTracking)
            sourceTracking = fullfile(dayDir,dirTracking(track_i).name);
            destTracking = [saveDir filesep dirBats(bat_i).name filesep dirBatDays(day_i).name filesep 'extracted' filesep];
            if ~exist(destTracking)
                mkdir(destTracking)
            end
            %copyfile(sourceTracking,destTracking);
        end
        
        %find analysis and processed dirs and copy
        dirBatSessions = dir('*_extraction');
        %dirBatSessions = dir([dirBats(bat_i).name(1:2) '*_extraction']);
        for session_i = 1:length(dirBatSessions)
            cd([dirBatSessions(session_i).name]);
            %find and copy over processed and analysis folders
            dirAnalysis = dir('analysis_*');
            dirProcessed = dir('processed_*');
            sourceAnalysis = fullfile(dirAnalysis(end).folder,dirAnalysis(end).name);
            sourceProcessed = fullfile(dirProcessed(end).folder,dirProcessed(end).name);
            destDir = [saveDir filesep dirBats(bat_i).name filesep dirBatDays(day_i).name filesep 'extracted' filesep dirBatSessions(session_i).name filesep];
            destAnalysis = fullfile(destDir,dirAnalysis(end).name);
            destProcessed = fullfile(destDir,dirProcessed(end).name);
            if ~exist(destAnalysis)
                mkdir(destAnalysis);
            end
            if ~exist(destProcessed)
                mkdir(destProcessed);
            end
            copyfile(sourceAnalysis,destAnalysis);
            copyfile(sourceProcessed,destProcessed);
            

            %go back to day directory
            cd(dayDir);
            disp(sourceAnalysis)
        end       
   cd(batDir);
   end
   cd(homeDir);
end