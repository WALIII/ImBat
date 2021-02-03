allBatsTag = 1;
batId = 'Ge';
saveFlag = 1;
saveTag = 'test';

if saveFlag == 1
    %saveDir1 = '/Volumes/Tobias_flig/topQualityData/analysis_done/plots/';
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\ForTobias\plots\';
    % Check if folder exists
    if exist([saveDir1 datestr(now,'yymmdd')])>0;
        disp('Youve been working today..');
    else
        mkdir([saveDir1 datestr(now,'yymmdd')])
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep];
end

if allBatsTag == 1
   allBats = [{'Gal'} {'Ge'} {'Ge_2'} {'Gi'} {'Go'} {'z2'} {'za'} {'zu'}];
   nBats = length(allBats);
else
    nBats = 1;
    
end

%initialize matrices
nDays = length(dayDir);
percPre = nan(nBats,nDays);
percPlace = nan(nBats,nDays);
percFlight = nan(nBats,nDays);
percPost = nan(nBats,nDays);
percAny = nan(nBats,nDays);
percNone = nan(nBats,nDays);


%for each day
for day_i = 1:length(dayDir)
   cd([dayDir(day_i).name filesep 'extracted']);
   %find all fly directories and enter last one
   flyDir = dir('*_fly-*_extraction');
   cd(flyDir(end).name);
   %find and enter last analysis folder
   analysisDir = dir('analysis_*');
   cd([analysisDir(end).name filesep 'placeCellsAng']);
   %load placeCell structure for % analysis
   placeCellStruct = dir('*_ExtractedPlaceCells_*');
   load(placeCellStruct.name);
   
   %save percentages into a cell
   percPre(day_i) = placeCells.perc_pre;
   percPlace(day_i) = placeCells.perc_place;
   percFlight(day_i) = placeCells.perc_place_loose;
   percPost(day_i) = placeCells.perc_post;
   percNone(day_i) = placeCells.perc_none;
   percAny(day_i) = 1 - placeCells.perc_none;
   
   %transfer the pre/place/post cells into respective folders
   ImBat_transferPlaceCells
   
   cd(dayDir(day_i).folder);
   
end
%concatenate and average all the days
percMean = [nanmean(percPre); nanmean(percPlace); nanmean(percPost); nanmean(percAny)];
percSEM = [nanstd(percPre,2)./sqrt(size(percPre,2)); nanstd(percPlace,2)./sqrt(size(percPlace,2)); nanstd(percPost,2)./sqrt(size(percPost,2)); nanstd(percAny,2)./sqrt(size(percAny,2))];
percAll = [percPre' percPlace' percFlight' percPost' percAny' percNone'];
%make figure for percentages
figPercMean = figure();
bar(1:nDays,percMean);
hold on;
errorbar(1:nDays,percMean,percSEM);
xlabel('Flight phase');
ylabel('%');
xticklabels({'Pre', 'Place', 'Post', 'Any'});
title([batId ': % of Pre, Place, Post Cells across ' num2str(nDays) ' days']);




