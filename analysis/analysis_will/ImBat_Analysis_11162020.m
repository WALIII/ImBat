function ImBat_Analysis_11162020(FlightAlignedROI)

% Cluster based on the difference of Flights
% prereq: [FlightAlignedROI] = ImBat_Align_FC(CombinedROI,flightPaths,2);

nclust = 3;
xyz2use = 4; % x = 1, y = 2; z = 3;


% behavioral analysis

VV_full = squeeze(FlightAlignedROI.ClustFlight(1:600,:,:));



if xyz2use ==1;
    t2use = ['Sorting based on X'];
    VV = squeeze(FlightAlignedROI.ClustFlight(1:600,xyz2use,:));
    
elseif xyz2use ==2
    t2use = ['Sorting based on Y'];
    VV = squeeze(FlightAlignedROI.ClustFlight(1:600,xyz2use,:));
    
elseif xyz2use ==3
    t2use = ['Sorting based on Z'];
    VV = squeeze(FlightAlignedROI.ClustFlight(1:600,xyz2use,:));
    
elseif xyz2use ==4
    t2use = ['Sorting based on PCA'];
    
    %% PCA!!
    VV_full2 = permute(VV_full,[1,3,2]); % permute order
    F2PCA = reshape(VV_full2,size(VV_full,1)*size(VV_full2,2),size(VV_full2,3)); % reshape matrix
    [coeff,score,latent,tsquared,explained,mu] = pca(F2PCA); % do PCA
    VV = reshape(score(:,1),size(VV_full,1),size(VV_full,3)); % reshape matrix
    
    explained(1); % report % explained...
    
end


CVC=corr(VV);

figure(); imagesc(CVC);
title( 'all-all flight consistancy');
xlabel('Flights');
ylabel('Flights');
colorbar
l = linkage(CVC, 'average', 'correlation');
c=cluster(l,'maxclust',nclust);
[aa,bb]=sort(c);


%% Color by time:
col1 = jet(size(VV_full,3));
titClust = {'X pos','Y pos','Z pos', '3D pos'};
% Plot the clusters
figure();
for ii = 1:3;
    subplot(1,4,ii)
    hold on;
    for i = 1:size(VV_full,3);
        plot((1: length(squeeze(VV_full(:,ii,i))))/120,squeeze(VV_full(:,ii,i)),'color',col1(i,:));
    end
    title(titClust(ii));
    xlabel('time (s)');
    ylabel('position, cm');
end
colormap(jet);
colorbar();
subplot(1,4,4)
title('Sorted by Time');

hold on;
for i = 1:size(VV_full,3);
    plot3(squeeze(VV_full(:,1,i)),squeeze(VV_full(:,2,i)),squeeze(VV_full(:,3,i)),'color',col1(i,:));
end


    col = jet(max(c));


if xyz2use <4;
    % resort c by the mean of the traj
    c_srt  = c;
    
    for i = 1: max(c);
        idx2use = find(c == i);
        %     plot3(squeeze(VV_full(:,1,idx2use)),squeeze(VV_full(:,2,idx2use)),squeeze(VV_full(:,3,idx2use)),'color',col(i,:));
        G(i) = mean(VV_full(340,xyz2use,idx2use),3);
    end
    
    c_srt = c;
    [a1 b1] = sort(G);
    for i = 1:max(c);
        idx2use = find(c == i);
        c_srt(idx2use) = b1(i);
    end
  CVC3=corr(VV(:,bb2));

%figure(); imagesc(corr(squeeze(FlightAlignedROI.ClustFlight(1:600,2,:)))); colorbar
figure(); imagesc(CVC3);
title( 'all-all flight consistancy');
xlabel('Flights');
ylabel('Flights');
colorbar


    
else
    c_srt = c;
end




%% Sort by XYZ or PC

    figure();
    
    % Plot the clusters
figure();
for ii = 1:3;
    subplot(1,4,ii)
    hold on;
    for i = 1: max(c);
        idx2use = find(c_srt == i);
        plot(squeeze(VV_full(:,ii,idx2use)),'color',col(i,:));
    end
end
title(t2use);
    subplot(1,4,4)

    hold on;
    for i = 1: max(c);
        idx2use = find(c_srt == i);
        plot3(squeeze(VV_full(:,1,idx2use)),squeeze(VV_full(:,2,idx2use)),squeeze(VV_full(:,3,idx2use)),'color',col(i,:));
        if i ==1;
            bb2 = idx2use;
        else
            bbtemp = idx2use;
            bb2 = cat(1, bb2,bbtemp);
            clear bbtemp
        end
    end



% max(c)
% col = jet(max(c));
% % Plot the clusters
% figure();
% hold on;
% for i = 1: max(c);
%     idx2use = find(c == i);
%     plot3(squeeze(VV_full(:,1,idx2use)),squeeze(VV_full(:,2,idx2use)),squeeze(VV_full(:,3,idx2use)),'color',col(i,:));
%
% end
% title(t2use);
%


PlotCells = 1;
if PlotCells ==1;
    figure();
    
    bound2plot = 1:500;
    
    CutCells = FlightAlignedROI.C_raw;
    disp(['Plotting ', num2str(size(CutCells,1)), ' Cells']);
    for i = 1:size(CutCells,1);%  for all cells
        cell2plot_unSorted = squeeze(CutCells(i,bound2plot,:))-mean(squeeze(CutCells(i,500:end,:)));
        cell2plot_Sorted = squeeze(CutCells(i,bound2plot,bb2))-mean(squeeze(CutCells(i,500:end,bb2)));
        subplot(1,2,1)
        imagesc(cell2plot_unSorted');
        title('sorted by days');
        subplot(1,2,2)
        
        imagesc(cell2plot_Sorted');
        title('sorted by cluster');
        
        pause();
        clf('reset');
    end
    
end



