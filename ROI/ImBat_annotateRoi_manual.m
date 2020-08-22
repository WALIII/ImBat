tifDir = dir('*_maxCorrProject.tif');
ROI = struct;
for i = 1:length(tifDir)
    ROI(i).batName = tifDir(i).name(1:3);
    ROI(i).dateSesh = tifDir(i).name(5:10);
    ROI(i).sessionType = tifDir(i).name(12:16);
    ROI(i).fileName = tifDir(i).name;
    ROI(i).folder = tifDir(i).folder;
    ROI(i).data = FS_annotate_image(tifDir(i).name);
end

for ii = length(tifDir):-1:1
        ROI(ii).data = FS_annotate_image(tifDir(ii).name);
end

save([ROI(1).batName '_ROI_' datestr(now,'YYmmDD_hhMM') '.mat'],'ROI');

%%
trackableROI = figure(); 
sgtitle([ROI(1).batName ': stable cells across ' num2str(length(ROI)) ' days']);
for iii = 1:length(ROI)
    subplot(4,3,iii);%,length(ROI.coordinates),iii); 
    hold on; 
    for p = 1:length(ROI(iii).data.coordinates) 
        plot(ROI(iii).data.coordinates{p}(:,1),ROI(iii).data.coordinates{p}(:,2)); 
        txt = text(ROI(iii).data.coordinates{p}(ceil(end/2),1),ROI(iii).data.coordinates{p}(1,2),num2str(p));
    end
    set(gca, 'YDir','reverse')
    title([ROI(iii).batName ' ' ROI(iii).dateSesh]);
    xticklabels([]);
    yticklabels([]);
    
end

savefig(trackableROI,['plot_' ROI(1).batName '_ROI_' datestr(now,'YYmmDD_hhMM') '.fig']);
saveas(trackableROI,['plot_' ROI(1).batName '_ROI_' datestr(now,'YYmmDD_hhMM') '.tif']);
