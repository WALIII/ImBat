function [out] = Deploy_SI(FlightAlignedROI)
% to calculate SI for a single day

% calculating Spatial Information
clear sigCell sigCell2 sigCell3 sigCell4

% params 
offset2use = 15; % frames of offset from flight to reward


for ii = 1:3
for i = 1:size(FlightAlignedROI{ii}.S,1); 
input_dat = squeeze(FlightAlignedROI{ii}.S(i,FlightAlignedROI{ii}.ROI_ON-10:FlightAlignedROI{ii}.ROI_ON+mean(FlightAlignedROI{ii}.FlightLength)/120*30+offset2use,:))';
input_dat(find(mean(input_dat')==0),:) = [];
 [sigCell(i,ii),~,~,tSI(i,ii)] = fTS_SI(input_dat);
 
  % pre and post pads
 input_dat2 = squeeze(FlightAlignedROI{ii}.S(i,FlightAlignedROI{ii}.ROI_ON-90:FlightAlignedROI{ii}.ROI_ON+mean(FlightAlignedROI{ii}.FlightLength)/120*30+90,:))';
input_dat2(find(mean(input_dat2')==0),:) = [];
 [sigCell2(i,ii),~,~,tSI2(i,ii)] = fTS_SI(input_dat2);

 % reward
 input_dat3 = squeeze(FlightAlignedROI{ii}.S(i,FlightAlignedROI{ii}.ROI_ON+mean(FlightAlignedROI{ii}.FlightLength)/120*30+offset2use:FlightAlignedROI{ii}.ROI_ON+mean(FlightAlignedROI{ii}.FlightLength)/120*30+90,:))';
input_dat3(find(mean(input_dat3')==0),:) = [];
 [sigCell3(i,ii),~,~,tSI3(i,ii)] = fTS_SI(input_dat3);
 
 % pre flight
  input_dat4 = squeeze(FlightAlignedROI{ii}.S(i,FlightAlignedROI{ii}.ROI_ON-90:FlightAlignedROI{ii}.ROI_ON-10,:))';
input_dat4(find(mean(input_dat4')==0),:) = [];
 [sigCell4(i,ii),~,~,tSI4(i,ii)] = fTS_SI(input_dat4);

 close all 
end
end

figure(); 
h2p = sum(sigCell');
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
out.Flight = sigCell;


clear h2p
figure(); 
h2p = sum(sigCell2');
h2p(h2p>0) = 1;
histogram(h2p,'Normalization','probability')
title('Including pre-post only');
ylim([0 1]);
 str={'','not-tuned','','Tuned',''};
 nTicks = numel( str );
xLim = get( gca, 'XLim' );
range = diff( xLim );
tickSpacing = range / nTicks;
xTicks = ( ( 1:nTicks  ) - 0.5 ) * tickSpacing + xLim(1);
set( gca, 'XTick', xTicks )
set( gca, 'XTickLabel', str )
out.All = sigCell2;

clear h2p
figure(); 
hold on;
h2p = sum(sigCell3');
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
out.Reward = sigCell3;


out.Pre = sigCell4;

