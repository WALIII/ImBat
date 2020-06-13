dirBat = dir('G*');
goodlist = strings(length(dirBat),5);
badlist = strings(length(dirBat),5);
badData = 0;

for i = 1:length(dirBat)
   cd(dirBat(i).name);
   dirSession = dir('*fly-*extraction');
%    goodlist(i) = char.empty;
%    badlist(i) = char.empty;
   for ii = 1:length(dirSession)
        cd(dirSession(ii).name);
        try
        dirAnalysis = dir('analysis*');
        analysisNewest = sort({dirAnalysis(end).name});
        analysisNewest = char(analysisNewest);
        %cd(analysisNewest);
        try
            openfig([dirAnalysis(1).folder filesep analysisNewest{1} filesep 'flights' filesep '*_flightVsVelocity.fig']);
        catch
            disp(['No cell data ' dirBat(i).name])
            badData = 1;
            badlist(i,ii) = dirSession(i).name;
        end
        try
            openfig([dirAnalysis(1).folder filesep analysisNewest{1} filesep 'ROI' filesep '*_maxProjectROI.fig']);
        catch
            disp(['No ROI data' dirBat(i).name])
            badData = 1;                
            badlist(i,ii) = dirSession(i).name;
        end
        if badData == o
            prompt = 'Is the imaging data good quality?';
            x = input(prompt)
            if x == 1
                goodlist(i,ii) = dirSession(i).name;
                %badlist(i,ii) = [];
            else
                badlist(i,ii) = dirSession(i).name;
                %goodlist(i,ii) = [];
            end
        end
        catch
        end
        close all;
        cd('../');
        badData = 0;
   end
   cd('../');
end