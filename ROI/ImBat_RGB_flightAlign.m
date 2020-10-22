function ImBat_RGB_flightAlign(batId,fullSeshTag,day1,day2,day3)
%batId = 'Gal';
%fullSeshTag = 1, chooses to look at whole session for 3 day comparison
clustNum = 2;
saveFlag = 0;
saveTag = 'clust2';
h = fspecial('disk',30);

%load first day placeCellStableROI data
if strcmp(batId,'Gal')
    %cd([dirTop(day1).folder filesep 'plots\200911-preDurPost cells across days']);
    load('Gal_200311to200324_activity_allTrials_allClusts_allTrials_sMat_newDff_newOrder_1to3.mat'); %load activity for pre,dur,post
elseif strcmp(batId,'Gen')
    load('Gen_200319to200324_activity_allTrials_allClusts_sMat_newDff.mat'); %load activity for pre,dur,post
end
%make saving directory
if saveFlag == 1
    %saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    saveDir1 = '/Users/periscope/Desktop/analysis/flight/plots/';
    if ~exist([saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign'])
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign']);
    else
        disp('You have been working today...');
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign' filesep];
end

if fullSeshTag == 0
    if isempty(day3) %comparing flight aligned max proj
        if day1==day2 %if looking at the cluster flight aligned vs full session
            %IM1 = mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day1});
            %IM2 = mat2gray(activity_allTrials.YmaxFull{day2});
            IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};
            IM2_raw = activity_allTrials.YmaxFull{day2};
            %normalize 2 images to each other
            IMcat = [IM1_raw IM2_raw]; %concatenate
            IMcat_gray = mat2gray(IMcat); %convert to 0-1
            IM1half = IMcat_gray(:,1:size(IMcat_gray,2)/2); %split into original matrices
            IM1preFilt = IM1half.*2; %multiply by 2 to equal the doubling of the cyan image
            IM1 = smoothdata(IM1preFilt,'gaussian',2);
            IM2preFilt = IMcat_gray(:,(size(IMcat_gray,2)/2)+1:end); %split into original matrices
            IM2 = smoothdata(IM2preFilt,'gaussian',2);
            plotTitle = [batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v full session (c)'];
        else %if looking at 1 day vs another day only
            %IM1 = mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day1});
            %IM2 = mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day2});
            IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};
            IM2_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2};
            %normalize the 2 images to each other
            IMcat = [IM1_raw IM2_raw];
            IMcat_gray = mat2gray(IMcat);
            IM1half = IMcat_gray(:,1:size(IMcat_gray,2)/2);
            IM1preFilt = IM1half.*2;
            IM1 = smoothdata(IM1preFilt,'gaussian',2);
            IM2preFilt = IMcat_gray(:,(size(IMcat_gray,2)/2)+1:end);
            IM2 = smoothdata(IM2preFilt,'gaussian',2);
            plotTitle = [batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (c)']
        end
        IM3 = [];
        % Alignment script
        rowTemp = [size(IM1,1),size(IM2,1),size(IM3,1)];
        colTemp = [size(IM1,2),size(IM2,2),size(IM3,2)];
        minRow = min(rowTemp(rowTemp>0));
        minCol = min(colTemp(colTemp>0));
        [IM1_aligned,IM2_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);
        IM3_aligned = [];
        % Will RGB script
        [a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);
        [a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);
    else %if looking at comparing 3 days
        IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};%mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day1});
        IM2_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2};%mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day2});
        IM3_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day3};%mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day3});
        IM1doub = imresize(double(IM1_raw),1);
        IM2doub = imresize(double(IM2_raw),1);
        IM3doub = imresize(double(IM3_raw),1);
        bground1=imfilter(IM1doub,h,'replicate');%smoothdata(IM1doub,'gaussian',2);%
        bground2=imfilter(IM2doub,h,'replicate');%smoothdata(IM2doub,'gaussian',2);%
        bground3=imfilter(IM3doub,h,'replicate');%smoothdata(IM3doub,'gaussian',2);%
        IM1_filt=IM1doub./(bground1+5);
        IM2_filt=IM2doub./(bground2+5);
        IM3_filt=IM3doub./(bground3+5);
        IMcat = [IM1_filt IM2_filt IM3_filt];
        IMcat_gray = mat2gray(IMcat);
        IM1 = IMcat_gray(:,1:size(IMcat_gray,2)/3);
        IM2 = IMcat_gray(:,(size(IMcat_gray,2)/3)+1:2*size(IMcat_gray,2)/3);
        IM3 = IMcat_gray(:,(2*size(IMcat_gray,2)/3)+1:end);
        
%         test=imresize(double(IM1preFilt),1);
%         x = [1 5 10 15 20 25];
%         % Plot
%         figure();
%         for i = 1: 6
%             subplot(2,3,i)
%             h=fspecial('disk',x(i));
%             bground=imfilter(test,h);
%             test=test./(bground+1);
%             imagesc(test);
%             title(['disk = ',num2str(x(i))]);
%         end

        %IM1 = smoothdata(IM1preFilt,'gaussian',2);
        %IM2 = smoothdata(IM2preFilt,'gaussian',2);
        %IM3 = smoothdata(IM3preFilt,'gaussian',2);
        %IM1 = imfilter(IM1preFilt,psf,'symmetric');
        %IM2 = imfilter(IM2preFilt,psf,'symmetric');
        %IM3 = imfilter(IM3preFilt,psf,'symmetric');
        minRow = min([size(IM1,1),size(IM2,1),size(IM3,1)]);
        minCol = min([size(IM1,2),size(IM2,2),size(IM3,2)]);
        % Alignment script
        [IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
        % Will RGB script
        [a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
        [a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);
        plotTitle = [batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
    end
elseif fullSeshTag == 1
    %if looking at comparing 3 days full sessions
    %IM1 = mat2gray(activity_allTrials.YmaxFull{day1});
    %IM2 = mat2gray(activity_allTrials.YmaxFull{day2});
    %IM3 = mat2gray(activity_allTrials.YmaxFull{day3});
    IM1_raw = activity_allTrials.YmaxFull{day1};
    IM2_raw = activity_allTrials.YmaxFull{day2};
    IM3_raw = activity_allTrials.YmaxFull{day3};
        IM1doub = imresize(double(IM1_raw),1);
        IM2doub = imresize(double(IM2_raw),1);
        IM3doub = imresize(double(IM3_raw),1);
        bground1=imfilter(IM1doub,h,'replicate');%smoothdata(IM1doub,'gaussian',2);%
        bground2=imfilter(IM2doub,h,'replicate');%smoothdata(IM2doub,'gaussian',2);%
        bground3=imfilter(IM3doub,h,'replicate');%smoothdata(IM3doub,'gaussian',2);%
        IM1_filt=IM1doub./(bground1+5);
        IM2_filt=IM2doub./(bground2+5);
        IM3_filt=IM3doub./(bground3+5);
        IMcat = [IM1_filt IM2_filt IM3_filt];    IMcat_gray = mat2gray(IMcat);
    IM1 = IMcat_gray(:,1:size(IMcat_gray,2)/3);
    IM2 = IMcat_gray(:,(size(IMcat_gray,2)/3)+1:2*size(IMcat_gray,2)/3);
    IM3 = IMcat_gray(:,(2*size(IMcat_gray,2)/3)+1:end);
    %IM1 = smoothdata(IM1preFilt,'gaussian',2);
    %IM2 = smoothdata(IM2preFilt,'gaussian',2);
    %IM3 = smoothdata(IM3preFilt,'gaussian',2);
    minRow = min([size(IM1,1),size(IM2,1),size(IM3,1)]);
    minCol = min([size(IM1,2),size(IM2,2),size(IM3,2)]);
    % Alignment script
    [IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
    % Will RGB script
    [a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
    [a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);
    plotTitle = [batId ' clust ' num2str(clustNum) ' Full Sesh: Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
end



plotDay123_unaligned = figure();
sgtitle(plotTitle);
subplot(1,2,1);
image((a1(:,:,:)));
title('un-aligned');
set(gca,'xticklabel',[],'yticklabel',[]);
subplot(1,2,2);
image((a2(:,:,:)));
title('aligned');
set(gca,'xticklabel',[],'yticklabel',[]);

plotDay123_aligned = figure();
image((a2(:,:,:)))
title(plotTitle);
set(gca,'xticklabel',[],'yticklabel',[]);
%set(gca,'YDir','normal');

if saveFlag == 1
    if strcmp(batId,'Gal')
        saveas(plotDay123_unaligned,[saveDir filesep 'Gal_200311and20_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_unaligned,[saveDir filesep 'Gal_200311and20_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_aligned,[saveDir filesep 'Gal_200311and20_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_aligned,[saveDir filesep 'Gal_200311and20_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
    elseif strcmp(batId,'Gen')
        saveas(plotDay123_unaligned,[saveDir filesep 'Gen_200319and24_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_unaligned,[saveDir filesep 'Gen_200319and24_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_aligned,[saveDir filesep 'Gen_200319and24_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_aligned,[saveDir filesep 'Gen_200319and24_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
    end
end

