function ImBat_PlotStoredDff(ROI_Data);

toPlot = size(ROI_Data,2);


figure();

for i = 1:toPlot
subplot(5,ceil(toPlot/5),i);
imagesc(ROI_Data{i}.MaxProj_flights);
title(num2str(ROI_Data{i}.date));
colormap(gray)
axis off
end


figure();

for i = 1:toPlot
subplot(5,ceil(toPlot/5),i);
imagesc(ROI_Data{i}.ROIs.results.Cn);
title(num2str(ROI_Data{i}.date));
colormap(gray)
axis off
end