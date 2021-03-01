function fTS_Aggregate_SI
% Deploy SI calculation across all bats:

% D02/18/21
% WAL3

% initialize vars
Reward1 = [];
All1 = [];
Flight1 = [];

DIR = pwd;
files = dir(pwd);
files(ismember( {files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.




    
    col = hsv(length(subFolders)+1);
    for i = 1 : length(subFolders)
        load([subFolders(i).name, '/CellReg_files/ROI_Data/Saved_Data/Aligned_Data.mat'],'FlightAlignedROI')
        %load([subFolders(i).name, '/CellReg_files/ROI_Data/ROI_Data.mat'])
        
        
        
        [out] = Deploy_SI(FlightAlignedROI);
        % combined ROIs
        Reward1 = cat(1,Reward1,out.Reward);
        Flight1 = cat(1,Flight1,out.Flight);
        All1 = cat(1,All1,out.All);
        % agg dat
        out_agg{i} = out;
        
        clear out

        
        clear FlightAlignedROI
        
    end
    
    
%% Plotting:
    figure(); 
h2p = sum(Flight1');
h2p(h2p>0) = 1;
histogram(h2p,'Normalization','probability')
title('Flight only');
ylim([0 1]);
 str={'','not-tuned','','Tuned',''};
 nTicks = numel( str );
xLim = get( gca, 'XLim' );
range = diff( xLim );
tickSpacing = range / nTicks;
xTicks = ( ( 1:nTicks  ) - 0.5 ) * tickSpacing + xLim(1);
set( gca, 'XTick', xTicks )
set( gca, 'XTickLabel', str )


clear h2p
figure(); 
h2p = sum(All1');
h2p(h2p>0) = 1;
histogram(h2p,'Normalization','probability')
title('Including pre-post ');
ylim([0 1]);
 str={'','not-tuned','','Tuned',''};
 nTicks = numel( str );
xLim = get( gca, 'XLim' );
range = diff( xLim );
tickSpacing = range / nTicks;
xTicks = ( ( 1:nTicks  ) - 0.5 ) * tickSpacing + xLim(1);
set( gca, 'XTick', xTicks )
set( gca, 'XTickLabel', str )


clear h2p
figure(); 
hold on;
h2p = sum(Reward1');
h2p(h2p>0) = 1;
histogram(h2p,'Normalization','probability')
title('Reward');
ylim([0 1]);

 str={'','not-tuned','','Tuned',''};
 nTicks = numel( str );
xLim = get( gca, 'XLim' );
range = diff( xLim );
tickSpacing = range / nTicks;
xTicks = ( ( 1:nTicks  ) - 0.5 ) * tickSpacing + xLim(1);
set( gca, 'XTick', xTicks )
set( gca, 'XTickLabel', str )

% plot bar graphs

for ii = 1: size(out_agg,2)
    
    h2p = sum(out_agg{ii}.Flight');
    h2p(h2p>0) = 1;
    FL(ii) = sum(h2p)/size(h2p,2); clear h2p
    
    h2p = sum(out_agg{ii}.All');
    h2p(h2p>0) = 1;
    AL(ii) = sum(h2p)/size(h2p,2); clear h2p
    
     h2p = sum(out_agg{ii}.Reward');
    h2p(h2p>0) = 1;
    RW(ii) = sum(h2p)/size(h2p,2); clear h2p
end

figure();
boxplot([FL' RW' AL'],'Labels',{'Flight','Reward','All'})
ylim([0 1]);   

  figure();
plotSpread([FL' RW' AL'])%,'Labels',{'Flight','Reward','All'})
ylim([0 1]);

