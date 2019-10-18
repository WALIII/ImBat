

function [CutCells, CutCells2,Idx] = ImBat_analysis_091219(ROI_Data);
% manually picked cells- how do things look all together?


bestCells{1} = [107, 58,60,68,34,95,9,29,47,41,16,26,132,67,98,148,20,94,137,115,122,75,111,92,127,83];
bestCells{2} = [16, 26, 107, 76, 74, 69, 122, 77,126,152,127,45,47,105,124,116,90,83,49];
bestCells{3} = [60, 75, 125, 64, 47,94, 29, 45, 69,122,131];
% build clustCells

for i = 1:3
[CutCells{i}, Ydata, ClustFlight{i},CutCells2{i}, CutFlights{i},Velocity{i}] = ImBat_Analysis_070119(ROI_Data,'day',12,'clust',i);
close all
end

figure(); hold on;  col = {'r','g','b'}; for i = 1:3 plot(CutFlights{i},col{i}); end
title( ' velocity profiles for clustered flights');
figure(); 
hold on;
plot3(squeeze(ClustFlight{1}(:,1,:)),squeeze(ClustFlight{1}(:,2,:)),squeeze(ClustFlight{1}(:,3,:)),'Color','r');
plot3(squeeze(ClustFlight{2}(:,1,:)),squeeze(ClustFlight{2}(:,2,:)),squeeze(ClustFlight{2}(:,3,:)),'Color','g');
plot3(squeeze(ClustFlight{3}(:,1,:)),squeeze(ClustFlight{3}(:,2,:)),squeeze(ClustFlight{3}(:,3,:)),'Color','b');
arrow3([squeeze(ClustFlight{1}(end-20,1,:)),squeeze(ClustFlight{1}(end-20,2,:)),squeeze(ClustFlight{1}(end-20,3,:))],[squeeze(ClustFlight{1}(end-10,1,:)),squeeze(ClustFlight{1}(end-10,2,:)),squeeze(ClustFlight{1}(end-10,3,:))]);
arrow3([squeeze(ClustFlight{2}(end-20,1,:)),squeeze(ClustFlight{2}(end-20,2,:)),squeeze(ClustFlight{2}(end-20,3,:))],[squeeze(ClustFlight{2}(end-10,1,:)),squeeze(ClustFlight{2}(end-10,2,:)),squeeze(ClustFlight{2}(end-10,3,:))]);
arrow3([squeeze(ClustFlight{3}(end-20,1,:)),squeeze(ClustFlight{3}(end-20,2,:)),squeeze(ClustFlight{3}(end-20,3,:))],[squeeze(ClustFlight{3}(end-10,1,:)),squeeze(ClustFlight{3}(end-10,2,:)),squeeze(ClustFlight{3}(end-10,3,:))]);

disp('oo');
%take 8 flights for each and concat the matrix for one in a subplot)
clear IM1 IM2 IM3

% subplot(1,3,1);
All_cells = cat(2, bestCells{1},bestCells{2},bestCells{3});

All_cells = unique(All_cells);


for Fl_sort = 1:3; 
%     Idx = bestCells{Fl_sort};
Idx = All_cells;
% Sort based on peak of mean for first cut...
[a b] = max(squeeze(mean(CutCells{Fl_sort}(Idx,50:300,:),3)),[],2);
[a1,b1] = sort(b);


for ii = 1:size(Idx,2)
    i = b1(ii);
    if ii ==1
IM1  =  zscore(squeeze(CutCells{1}(Idx(i),:,1:9)));
IM2  =  zscore(squeeze(CutCells{2}(Idx(i),:,1:9)));
IM3  =  zscore(squeeze(CutCells{3}(Idx(i),:,1:9)));
% mean cells
for iii = 1:3
mIM(:,:,iii) = zscore(squeeze(mean(CutCells{iii}(Idx(b1),1:578,:),3)),[],2);
end
    else
        
IM1  = cat(2,IM1,zeros(size(IM1,1),1)+4);
IM2  = cat(2,IM2,zeros(size(IM2,1),1)+4);
IM3  = cat(2,IM3,zeros(size(IM3,1),1)+4);

IM1  = cat(2,IM1,zscore(squeeze(CutCells{1}(Idx(i),:,1:9))));
IM2  = cat(2,IM2,zscore(squeeze(CutCells{2}(Idx(i),:,1:9))));
IM3  = cat(2,IM3,zscore(squeeze(CutCells{3}(Idx(i),:,1:9))));
    end
end
figure();
subplot(131)
imagesc(IM1',([0 4])); 
subplot(132)
imagesc(IM2',([0 4])); 
subplot(133)
imagesc(IM3',([0 4])); 
colormap(hot)
figure(); [RGB1 RGB2] = CaBMI_XMASS(squeeze(mIM(:,1:578,1)),squeeze(mIM(:,1:578,2)),squeeze(mIM(:,1:578,2)),'HL',[0.3 0.7]); figure(); imagesc(squeeze(RGB1(:,:,1,:)))

figure(); [RGB1 RGB2] = CaBMI_XMASS(IM1(1:578,:)',IM2(1:578,:)',IM3(1:578,:)','HL',[0.3 0.7]); figure(); imagesc(squeeze(RGB1(:,:,1,:)))
end


