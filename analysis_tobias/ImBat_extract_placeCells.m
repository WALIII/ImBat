allBatsTag = 0;
inspectNoTuneFlag = 0;
batId = 'Ge';
saveFlag = 1;
saveTag = 'mean-bitsPerSec-max-noSignifPlace';

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

nDays = 2;%length(dayDir);
%initialize matrices
percPre = nan(nBats,nDays);
percPlace = nan(nBats,nDays);
percFlight = nan(nBats,nDays);
percPost = nan(nBats,nDays);
percAny = nan(nBats,nDays);
percNone = nan(nBats,nDays);



%for each day, calculate percentages of each phase
for day_i = 1:2%length(dayDir)
   cd([dayDir(day_i).name filesep 'extracted']);
   %find all fly directories and enter last one
   flyDir = dir('*_fly-*_extraction');
   cd(flyDir(end).name);
   %find and enter last analysis folder
   analysisDir = dir('analysis_*');
   cd([analysisDir(end-4).name filesep 'placeCellsAng']);
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
   [placeCells ] = ImBat_transferPlaceCells(placeCells)
  
   %find all the cells not tuned to any phase
   [placeCells] = ImBat_extract_noResponseCells(placeCells)
   cd ..;
   save(placeCellStruct.name,'placeCells');
   cd(dayDir(day_i).folder);
   
   if inspectNoTuneFlag == 1
       % ImBat_inspect_noTuneCells(placeCells)
   end
   
end

%% concatenate and average all the days
percMean = [nanmean(percPre); nanmean(percPlace); nanmean(percPost); nanmean(percAny)];
percSEM = [nanstd(percPre)./sqrt(size(percPre,2)); nanstd(percPlace)./sqrt(size(percPlace,2)); nanstd(percPost)./sqrt(size(percPost,2)); nanstd(percAny)./sqrt(size(percAny,2))];
percAll = [percPre' percPlace' percFlight' percPost' percAny' percNone'];
%make figure for percentages
figPercMean = figure();
bar(percMean);
hold on;
er = errorbar(percMean,percSEM);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xlabel('Flight phase');
ylabel('%');
xticklabels({'Pre', 'Place', 'Post', 'Any'});
title([batId ': Mean % of Pre, Place, Post Cells across ' num2str(nDays) ' days ' saveTag]);

figPercAll = figure();
bar(percAll);
hold on;
xlabel('Day');
ylabel('%');
legend({'Pre', 'Place','Flight','Post', 'Any','None'});
title([batId ': % of Pre, Place, Post Cells across ' num2str(nDays) ' days ' saveTag]);

%save figs
if saveFlag == 1
savefig(figPercMean,[saveDir batId '_percPrePlacePostCells_mean_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
saveas(figPercMean, [saveDir batId '_percPrePlacePostCells_mean_' saveTag '_' datestr(now,'yymmdd-HHMM') '.jpg']);
savefig(figPercAll,[saveDir batId '_percPrePlacePostCells_all_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
saveas(figPercAll, [saveDir batId '_percPrePlacePostCells_all_' saveTag '_' datestr(now,'yymmdd-HHMM') '.jpg']);
end




