batId = 'Gen';
saveFlag = 1;
nFlights = 
clustNum = 2;
lenTrace = length(dataPreDurPost.mean_act_aligned{clustNum}{1,1});
nDays = length(dataPreDurPost.mean_act_aligned{clustNum}(:,1));
%make saving directory
if saveFlag == 1
    saveDir1 = '/Volumes/Tobias_flig/topQualityData/analysis_done/plots/';
    if ~exist([saveDir1 datestr(now,'yymmdd') filesep 'flightAlignMaxProj'])
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'flightAlignMaxProj']);
    else
        disp('You have been working today...');
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'flightAlignMaxProj' filesep];
end

framesFlightAlign = zeros(m,n,nFlights)