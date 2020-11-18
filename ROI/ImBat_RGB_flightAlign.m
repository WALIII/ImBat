function ImBat_RGB_flightAlign(batId,fullSeshTag,day1,day2,day3)
%batId = 'Gal';
%fullSeshTag = 1, chooses to look at whole session for 3 day comparison
clustNum = 2;
saveFlag = 1;
rigidFlag = 1;
saveTag = 'rigid';
h = fspecial('disk',50);

%load first day placeCellStableROI data
if strcmp(batId,'Gal')
    %cd([dirTop(day1).folder filesep 'plots\200911-preDurPost cells across days']);
    load('Gal_200311to200324_activity_allTrials_allClusts_allTrials_sMat_newDff_newOrder.mat'); %load activity for pre,dur,post
elseif strcmp(batId,'Gen')
    load('Gen_200319to200324_activity_allTrials_allClusts_sMat_newDff_newOrder.mat'); %load activity for pre,dur,post
end
%make saving directory
if saveFlag == 1
    %saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    %saveDir1 = '/Users/periscope/Desktop/analysis/flight/plots/';
    saveDir1 = 'C:\Users\tobias\Desktop\analysis\plots\';
    if ~exist([saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign'])
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign']);
    else
        disp('You have been working today...');
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign' filesep];
end

% %set raw to the flight aligned frames,
% IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};
% IM2_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2};
% IM1_rawFull = activity_allTrials.YmaxFull{day1};
% IM2_rawFull = activity_allTrials.YmaxFull{day2};
% if ~isempty(day3) %if there is a day3
%     IM3_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day3};
%     IM3_rawFull = activity_allTrials.YmaxFull{day3};
% else
%     IM3_raw = IM2_raw; %replicate for now
%     IM3_rawFull = IM2_rawFull;
% end
%pull out the correct data into the raw image variables
if fullSeshTag == 0 %comparing flight aligned max proj
    if isempty(day3)
        if day1 == day2 %if looking at the cluster flight aligned vs full session
            %set raw to the flight aligned frames and full session
            IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};
            IM2_raw = activity_allTrials.YmaxFull{day1};
            IM3_raw = activity_allTrials.YmaxFull{day1}; %duplicate for now
            plotTitle = [saveTag ' ' batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v full session (c)'];
        else %if looking at 1 day vs another day only
            %set raw to the flight aligned frames,
            IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};
            IM2_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2};
            IM3_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2}; %duplicate for now
            plotTitle = [saveTag ' ' batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (c)']
        end
    else %if looking at comparing 3 days
        IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};
        IM2_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2};
        IM3_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day3};
        plotTitle = [saveTag ' ' batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
     end
else %comparing the full session
    if isempty(day3) %comparing 1 day vs another
        IM1_raw = activity_allTrials.YmaxFull{day1};
        IM2_raw = activity_allTrials.YmaxFull{day2};
        IM3_raw = activity_allTrials.YmaxFull{day2}; %duplicate for now
        plotTitle = [saveTag ' ' batId ' clust ' num2str(clustNum) ' Full Sesh: Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (c)'];
    else %comparing 3 days full
        IM1_raw = activity_allTrials.YmaxFull{day1};
        IM2_raw = activity_allTrials.YmaxFull{day2};
        IM3_raw = activity_allTrials.YmaxFull{day3};
        plotTitle = [saveTag ' ' batId ' clust ' num2str(clustNum) ' Full Sesh: Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
    end
end


%filter, and divide background to make the pixels all the same brightness
%resize and subtract background
IM1doub = imresize(double(IM1_raw),2);
IM2doub = imresize(double(IM2_raw),2);
IM3doub = imresize(double(IM3_raw),2);
bground1=imfilter(IM1doub,h,'replicate');%smoothdata(IM1doub,'gaussian',2);%
bground2=imfilter(IM2doub,h,'replicate');%smoothdata(IM2doub,'gaussian',2);%
bground3=imfilter(IM3doub,h,'replicate');%smoothdata(IM3doub,'gaussian',2);%
IM1_filt=IM1doub./(bground1+5);
IM2_filt=IM2doub./(bground2+5);
IM3_filt=IM3doub./(bground3+5);
%concatenate to make into grayscale and then break apart into correct sizes
IMcat = [IM1_filt IM2_filt IM3_filt];    IMcat_gray = mat2gray(IMcat);
IM1 = IMcat_gray(:,1:size(IMcat_gray,2)/3);
IM2 = IMcat_gray(:,(size(IMcat_gray,2)/3)+1:2*size(IMcat_gray,2)/3);
IM3 = IMcat_gray(:,(2*size(IMcat_gray,2)/3)+1:end);
%make sure all images are same size if overlapping later
minRow = min([size(IM1,1),size(IM2,1),size(IM3,1)]);
minCol = min([size(IM1,2),size(IM2,2),size(IM3,2)]);
%make copy of the image to heavily filter for image registration
image1 = IM1;
image2 = IM2;
image3 = IM3;
% Clip image ( binarize) above 25 percentile
image1(image1<0.25) =0;
image2(image2<0.25) =0;
image3(image3<0.25) =0;
% Clip image ( binarize) above 25 percentile
image1(image1>=0.75) =1;
image2(image2>=0.75) =1;
image3(image3>=0.75) =1;
% remove edges ( in case these are producing artifacts)
edgeCutoff = 60;
image1(1:edgeCutoff,:) = 0;
image1(:,1:edgeCutoff) = 0;
image1(end-edgeCutoff:end,:) = 0;
image1(:,end-edgeCutoff:end) = 0;
image2(1:edgeCutoff,:) = 0;
image2(:,1:edgeCutoff) = 0;
image2(end-edgeCutoff:end,:) = 0;
image2(:,end-edgeCutoff:end) = 0;
image3(1:edgeCutoff,:) = 0;
image3(:,1:edgeCutoff) = 0;
image3(end-edgeCutoff:end,:) = 0;
image3(:,end-edgeCutoff:end) = 0;
if isempty(day3)
    IM3 = [];
    image3 = [];
end
if rigidFlag == 1 %use the rigid alignment with affine
    try
        [IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(image1,image2,image3,IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
    catch
        [IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(image1,image2,image3,IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);
    end
    [a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);
else %use nonrigid with demon flow
    [figNonrigid,imagesAligned] = ImBat_imageAlign_nonrigid(image1,image2,image3,IM1,IM2,IM3);
    [a2,b2] = CaBMI_XMASS(imagesAligned.IM1_aligned,imagesAligned.IM2_aligned,imagesAligned.IM3_aligned);
end
[a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);

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
            if rigidFlag == 0
                saveas(figNonrigid,[saveDir filesep 'Gal_200311and20_nonrigidAlign_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
                savefig(figNonrigid,[saveDir filesep 'Gal_200311and20_nonrigidAlign_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
            end
        elseif strcmp(batId,'Gen')
            saveas(plotDay123_unaligned,[saveDir filesep 'Gen_200319and24_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
            savefig(plotDay123_unaligned,[saveDir filesep 'Gen_200319and24_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
            saveas(plotDay123_aligned,[saveDir filesep 'Gen_200319and24_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
            savefig(plotDay123_aligned,[saveDir filesep 'Gen_200319and24_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
            if rigidFlag == 0
                saveas(figNonrigid,[saveDir filesep 'Gen_200319and24_nonrigidAlign_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
                savefig(figNonrigid,[saveDir filesep 'Gen_200319and24_nonrigidAlign_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
            end
        end
    end
    
    
%%

%     
%     if fullSeshTag == 0
%         if isempty(day3) %comparing flight aligned max proj
%             if day1==day2 %if looking at the cluster flight aligned vs full session
%                 
%                 plotTitle = [batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v full session (c)'];
%             else %if looking at 1 day vs another day only
%                 IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};
%                 IM2_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2};
%                 %normalize the 2 images to each other
%                 IMcat = [IM1_raw IM2_raw];
%                 IMcat_gray = mat2gray(IMcat);
%                 IM1half = IMcat_gray(:,1:size(IMcat_gray,2)/2);
%                 IM1preFilt = IM1half.*2;
%                 IM1 = smoothdata(IM1preFilt,'gaussian',2);
%                 IM2preFilt = IMcat_gray(:,(size(IMcat_gray,2)/2)+1:end);
%                 IM2 = smoothdata(IM2preFilt,'gaussian',2);
%                 plotTitle = [batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (c)']
%             end
%             IM3 = [];
%             % Alignment script
%             rowTemp = [size(IM1,1),size(IM2,1),size(IM3,1)];
%             colTemp = [size(IM1,2),size(IM2,2),size(IM3,2)];
%             minRow = min(rowTemp(rowTemp>0));
%             minCol = min(colTemp(colTemp>0));
%             [IM1_aligned,IM2_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);
%             IM3_aligned = [];
%             % Will RGB script
%             [a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);
%             [a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);
%         else %if looking at comparing 3 days
%             IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};%mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day1});
%             IM2_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2};%mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day2});
%             IM3_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day3};%mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day3});
%             %double size of data so everything isn't washed out by removing background
%             IM1doub = imresize(double(IM1_raw),2);
%             IM2doub = imresize(double(IM2_raw),2);
%             IM3doub = imresize(double(IM3_raw),2);
%             %find background of the image with a fspecial filter from above
%             bground1=imfilter(IM1doub,h,'replicate');%smoothdata(IM1doub,'gaussian',2);%
%             bground2=imfilter(IM2doub,h,'replicate');%smoothdata(IM2doub,'gaussian',2);%
%             bground3=imfilter(IM3doub,h,'replicate');%smoothdata(IM3doub,'gaussian',2);%
%             %remove the background and amplify low intensity pixels so they are same power as high intensity pixels
%             IM1_filt=IM1doub./(bground1+1);
%             IM2_filt=IM2doub./(bground2+1);
%             IM3_filt=IM3doub./(bground3+1);
%             %convert to grayscale and cut back apart into correct segments
%             IMcat = [IM1_filt IM2_filt IM3_filt];
%             IMcat_gray = mat2gray(IMcat);
%             IM1 = IMcat_gray(:,1:size(IMcat_gray,2)/3);
%             IM2 = IMcat_gray(:,(size(IMcat_gray,2)/3)+1:2*size(IMcat_gray,2)/3);
%             IM3 = IMcat_gray(:,(2*size(IMcat_gray,2)/3)+1:end);
%             %         test=imresize(double(IM1preFilt),1);
%             %         x = [1 5 10 15 20 25];
%             %         % Plot
%             %         figure();
%             %         for i = 1: 6
%             %             subplot(2,3,i)
%             %             h=fspecial('disk',x(i));
%             %             bground=imfilter(test,h);
%             %             test=test./(bground+1);
%             %             imagesc(test);
%             %             title(['disk = ',num2str(x(i))]);
%             %         end
%             minRow = min([size(IM1,1),size(IM2,1),size(IM3,1)]);
%             minCol = min([size(IM1,2),size(IM2,2),size(IM3,2)]);
%             % Alignment script
%             [IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
%             % Will RGB script
%             [a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
%             [a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);
%             plotTitle = [batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
%         end
%     elseif fullSeshTag == 1
%         %if looking at comparing 3 days full sessions
%         %IM1 = mat2gray(activity_allTrials.YmaxFull{day1});
%         %IM2 = mat2gray(activity_allTrials.YmaxFull{day2});
%         %IM3 = mat2gray(activity_allTrials.YmaxFull{day3});
%         IM1_raw = activity_allTrials.YmaxFull{day1};
%         IM2_raw = activity_allTrials.YmaxFull{day2};
%         IM3_raw = activity_allTrials.YmaxFull{day3};
%         %IM1_scale = imresize(IM1_raw,1);
%         %IM2_scale = imresize(IM2_raw,1);
%         %IM3_scale = imresize(IM3_raw,1);
%         IM1doub = imresize(double(IM1_raw),2);
%         IM2doub = imresize(double(IM2_raw),2);
%         IM3doub = imresize(double(IM3_raw),2);
%         bground1=imfilter(IM1doub,h,'replicate');%smoothdata(IM1doub,'gaussian',2);%
%         bground2=imfilter(IM2doub,h,'replicate');%smoothdata(IM2doub,'gaussian',2);%
%         bground3=imfilter(IM3doub,h,'replicate');%smoothdata(IM3doub,'gaussian',2);%
%         IM1_filt=IM1doub./(bground1+5);
%         IM2_filt=IM2doub./(bground2+5);
%         IM3_filt=IM3doub./(bground3+5);
%         IMcat = [IM1_filt IM2_filt IM3_filt];    IMcat_gray = mat2gray(IMcat);
%         IM1 = IMcat_gray(:,1:size(IMcat_gray,2)/3);
%         IM2 = IMcat_gray(:,(size(IMcat_gray,2)/3)+1:2*size(IMcat_gray,2)/3);
%         IM3 = IMcat_gray(:,(2*size(IMcat_gray,2)/3)+1:end);
%         %IM1 = smoothdata(IM1preFilt,'gaussian',2);
%         %IM2 = smoothdata(IM2preFilt,'gaussian',2);
%         %IM3 = smoothdata(IM3preFilt,'gaussian',2);
%         minRow = min([size(IM1,1),size(IM2,1),size(IM3,1)]);
%         minCol = min([size(IM1,2),size(IM2,2),size(IM3,2)]);
%         image1 = IM1;
%         image2 = IM2;
%         image3 = IM3;
%         % Clip image ( binarize) above 25 percentile
%         image1(image1<0.25) =0;
%         image2(image2<0.25) =0;
%         image3(image3<0.25) =0;
%         % Clip image ( binarize) above 25 percentile
%         image1(image1>=0.75) =1;
%         image2(image2>=0.75) =1;
%         image3(image3>=0.75) =1;
%         % remove edges ( in case these are producing artifacts)
%         image1(1:10,1:10) = 0;
%         image2(1:10,1:10) = 0;
%         image3(1:10,1:10) = 0;
%         image1(end-10:end,end-10:end) = 0;
%         image2(end-10:end,end-10:end) = 0;
%         image3(end-10:end,end-10:end) = 0;
%         
%         [figNonrigid,imagesAligned] = ImBat_imageAlign_nonrigid(image1,image2,image3,IM1,IM2,IM3);
%         % Alignment script
%         %[IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
%         % Will RGB script
%         [a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
%         [a2,b2] = CaBMI_XMASS(imagesAligned.IM1_aligned,imagesAligned.IM2imagesAligned.,IM3_aligned);
%         plotTitle = ['nonrigid alpha' num2str(imagesAligned.alpha) ' ' batId ' clust ' num2str(clustNum) ' Full Sesh: Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
%     end
%     
%     
%     
%    
    
