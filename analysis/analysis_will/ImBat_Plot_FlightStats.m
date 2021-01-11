function ImBat_Plot_FlightStats(out,col);

% start with: out = ImBat_FlightStats(FlightAlignedROI);
i = length(out.toPlot);

hold on;
 plot(out.toPlot(1,:),'color',col,'LineWidth',2);
sz = out.toPlot(3,:);
sz(sz ==0) = nan;
s = scatter(1:length(out.toPlot(1,:)),out.toPlot(1,:),sz*5,col,'filled');
s.MarkerFaceAlpha = 0.4;
    y1 = out.toPlot(1,:);
    x1 = 1:length(y1);
    err1 = out.toPlot(2,:);
h = errorbar(x1,y1,err1,'vertical','color',col,'LineWidth',1.2);%,'CapSize',1);

% Set transparency level (0:1)
alpha = 0.4;   
% Set transparency (undocumented)
set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
% figure(); plot(pval_combined_data);

% [a b] = find(out.pval_combined_data<(0.05./length(out.toPlot)));
% scatter(b,ones(length(b),1),'k*')
title('flight correlation vs time');
ylabel(' Correlation to day 1');
xlabel(' days');
