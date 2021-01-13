dirGbats = dir('G*');
dirZbats = dir('z*');
saveDir = 'F:\flight\topQualityData\copydir';
for bat_i = 1:length(dirGbats)
   cd(dirGbats(bat_i).name);
   dirBatDays = dir([dirGbats(bat_i).name(1:2) '*']);
   for day_i = 1:length(dirBatDays)
        cd([dirBatDays(day_i).name filesep 'extracted']);
        dirBatSessions = dir([dirGbats(bat_i).name(1:2) '*_extraction']);
        for session_i = 1:length(dirBatSessions)
            cd([dirBatSessions(session_i).name]);
            dirAnalysis = dir('analysis_*');
            dirProcessed = dir('processed_*');
            sourceAnalysis = fullfile(dirAnalysis(end).folder,dirAnalysis(end).name);
            sourceProcessed = fullfile(dirProcessed(end).folder,dirProcessed(end).name);
            destDir = [saveDir filesep dirGbats(bat_i).name filesep dirBatDays(day_i).name filesep 'extracted' filesep dirBatSessions(session_i).name filesep];
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
        end
        
   
   end
end