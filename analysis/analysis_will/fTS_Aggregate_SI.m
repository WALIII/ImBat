function fTS_Aggregate_SI
% Deploy SI calculation across all bats:

% D02/18/21
% WAL3

warning off

% initialize vars
Reward1 = [];
All1 = [];
Flight1 = [];
Pre1 = [];
FL_total(1:5,:) = 0;


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
        Pre1 = cat(1,Pre1,out.Pre);

        
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
        FL_total(1,:) = FL_total(1,:)+sum(h2p); % FL_total has all
        FL_total(5,:) = FL_total(5,:)+size(h2p,2); % SIZE

    FL(ii) = sum(h2p)/size(h2p,2); clear h2p
    
    h2p = sum(out_agg{ii}.All');
    h2p(h2p>0) = 1;
        FL_total(2,:) = FL_total(2,:)+sum(h2p);
    AL(ii) = sum(h2p)/size(h2p,2); clear h2p
    
     h2p = sum(out_agg{ii}.Reward');
    h2p(h2p>0) = 1;
        FL_total(3,:) = FL_total(3,:)+sum(h2p);
    RW(ii) = sum(h2p)/size(h2p,2); clear h2p
    
    h2p = sum(out_agg{ii}.Pre');
    h2p(h2p>0) = 1;
        FL_total(4,:) = FL_total(4,:)+sum(h2p);
    PR(ii) = sum(h2p)/size(h2p,2); clear h2p
end

figure();
boxplot([PR' FL' RW' AL'],'Labels',{'Pre' 'Flight','Reward','All'})
ylim([0 1]);   

  figure();
plotSpread([PR' FL' RW' AL'],'showMM',4)%,'Labels',{'Flight','Reward','All'})
ylim([0 1]);

FL_total2 = round(FL_total*0.7087); % flight, all_flight_phase, reward, pre, all_cells

