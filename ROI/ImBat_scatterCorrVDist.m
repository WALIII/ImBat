figure(); 
for i = 1:4
   subplot(2,2,i)
    scatter(ROI_refined.timeCorrS{i}(ROI_refined.corrIndAllS{i}),ROI_refined.distance{i}(ROI_refined.corrIndAll{i})) 
    hold on;
    sgtitle('Smatrix: Pairwise ROI distance (um) vs corr coeff (R)')
    title(['Gal Day# ' num2str(i)]);
    if i == 3
        xlabel('CorrCoef');
        ylabel('Distance (um)');
    end
end

figure(); 
for i = 1:4
   subplot(2,2,i)
    scatter(ROI_refined.timeCorrCloseS{i}(ROI_refined.corrIndCloseS{i}),ROI_refined.distance{i}(ROI_refined.corrIndClose{i}),'r') 
    hold on;
    sgtitle('SMatrix: Pairwise close ROI distance (um) vs corr coeff (R)')
    title(['Gal Day# ' num2str(i)]);
    if i == 3
        xlabel('CorrCoef');
        ylabel('Distance (um)');
    end
end