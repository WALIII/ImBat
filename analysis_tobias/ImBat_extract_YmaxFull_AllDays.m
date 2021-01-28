function [YmaxFullAllDays] = ImBat_extract_YmaxFull_AllDays(batId)
saveFlag = 1; %do you want to save the figures and output structure?
rest1Flag = 1; %do you want to include the rest1 session too?

saveTag = 'rest';
if saveFlag == 1
    %saveDir1 = '/Volumes/Tobias_flig/topQualityData/analysis_done/plots/';
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\ForTobias\plots\';
    % Check if folder exists
    if exist([saveDir1 datestr(now,'yymmdd')])>0;
        disp('Youve been working today..');
    else
        mkdir([saveDir1 datestr(now,'yymmdd')])
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep];
end

if strcmp(batId,'Gal')
    nDays = [1:40]; %which days to look at
    dirTop = dir('Ga*');
elseif strcmp(batId,'Gen')
    nDays = [1:100];
    dirTop = dir('Ge*');
elseif strcmp(batId,'Z2')
    nDays = [1:40]; %which days to look at
    dirTop = dir('z2*');
elseif strcmp(batId, 'Zu')
    nDays = [1:40];
    dirTop = dir('zu*');
elseif strcmp(batId,'Za')
    nDays = [1:5];
    dirTop = dir('za*');
end


for day_i = 1:length(nDays) %[50,51,66,67,76,77]%[2,3,27,28,33,34]%for each day
    %load results data
    try %extract metadata names and enter processed folder
        cd([dirTop(nDays(day_i)).name filesep 'extracted'])
        flyFolders = dir('*fly*extraction');
        rest1Folders = dir('*rest1*extraction');
        batName{day_i} = flyFolders(end).name(1:4);
        dateSesh{day_i} = flyFolders(end).name(6:11);
        sessionType{day_i} = flyFolders(end).name(13:17);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
            cd(dirProcessed(end).name); %can change this if need to look at earlier or later processed folders based on batname, date, etc
    catch
        cd(dirTop(nDays(day_i)).name);
        flyFolders = dir('*fly*extraction');
        rest1Folders = dir('*rest1*extraction');
        batName{day_i} = flyFolders(end).name(1:3);%(1:3);
        dateSesh{day_i} = flyFolders(end).name(5:10);%(5:10);
        sessionType{day_i} = flyFolders(end).name(12:16); %(12:16);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        cd(dirProcessed(end).name);
    end
    %load('results.mat'); %load cellData
    vidData = load('Motion_corrected_Data_DS.mat'); %load video frame data
    %make gaussian filter based on dff function for spatial filtering
    gSig = 1;
    gSiz = 4.5*gSig;
    scaling = 10; %scale up the ymax to make less pixelated
    psf = fspecial('gaussian', round(2*gSiz), gSig);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;   % only use pixels within the center disk
    %   % Take median, filter, and max of full movie
    Y_med = median(vidData.Y,3);
    Y_min = min(vidData.Y,3);
    Ydff = vidData.Y - Y_med; %subtract med
    Ydff_tFilt = medfilt3(Ydff); %temporal filtering
    Ydff_filt = imfilter(Ydff_tFilt,psf,'symmetric'); %spatial filtering
    length_vector = 1:length(Ydff_filt);
    rand_length_vector= randperm(length(length_vector));
    YmaxFull{15} = max(Ydff_filt(:,:,rand_length_vector(1:round(length(length_vector)/3))),[],3);
    YmaxFull{16} = max(Ydff_filt(:,:,rand_length_vector(round(length(length_vector)/3):round(length(length_vector)/3)*2)),[],3);
    YmaxFull{17} = max(Ydff_filt(:,:,rand_length_vector(round(length(length_vector)/3)*2:end)),[],3);
    
    YmaxFull{day_i} = max(Ydff_filt,[],3); %take max
    YmaxFull{day_i} = imresize(YmaxFull{day_i},scaling); %resize to eliminate pixelation
    %alignment = load('Alignment.mat');
    cd(dirProcessed(end).folder);
    disp(num2str(day_i))
    clear vidData;
    cd ..
    if rest1Flag == 1
        cd(rest1Folders(end).name);
        dirProcessed = dir('processed_*');
            cd(dirProcessed(end).name); %can change this if need to look at earlier or later processed folders based on batname, date, etc
       %load('results.mat'); %load cellData
    vidData = load('Motion_corrected_Data_DS.mat'); %load video frame data
    %   % Take median, filter, and max of full movie
    Y_med = median(vidData.Y,3);
    Y_min = min(vidData.Y,3);
    Ydff = vidData.Y - Y_med; %subtract med
    Ydff_tFilt = medfilt3(Ydff); %temporal filtering
    Ydff_filt = imfilter(Ydff_tFilt,psf,'symmetric'); %spatial filtering
    length_vector = 1:length(Ydff_filt);
    rand_length_vector= randperm(length(length_vector));
    YmaxRest{15} = max(Ydff_filt(:,:,rand_length_vector(1:round(length(length_vector)/3))),[],3);
    YmaxRest{16} = max(Ydff_filt(:,:,rand_length_vector(round(length(length_vector)/3):round(length(length_vector)/3)*2)),[],3);
    YmaxRest{17} = max(Ydff_filt(:,:,rand_length_vector(round(length(length_vector)/3)*2:end)),[],3);
    
    YmaxRest{day_i} = max(Ydff_filt,[],3); %take max
    YmaxRest{day_i} = imresize(YmaxRest{day_i},scaling); %resize to eliminate pixelation
    %alignment = load('Alignment.mat');
    cd(dirProcessed(end).folder);
    disp(num2str(day_i))
    clear vidData;
    
    
    end
    
    cd(dirTop(1).folder);

end

YmaxFullAllDays.batId = batId;
YmaxFullAllDays.batName = batName;
YmaxFullAllDays.sessionType = sessionType;
YmaxFullAllDays.dateSesh = dateSesh;
YmaxFullAllDays.YmaxFull = YmaxFull;
if rest1Flag == 1
    YmaxFullAllDays.YmaxRest = YmaxRest;
end

if saveFlag == 1
    if strcmp(batId,'Gal')
        save([saveDir 'Gal_200227to200404_YmaxFull_' saveTag '.mat'],'YmaxFullAllDays','-v7.3');
    elseif strcmp(batId,'Gen')
        save([saveDir 'Gen_200305to200311_YmaxFull_' saveTag '.mat'],'YmaxFullAllDays','-v7.3');
        %save([saveDir 'Gen_200424to200523_YmaxFull_' saveTag '.mat'],'YmaxFullAllDays','-v7.3');
    elseif strcmp(batId,'Z2')
        save([saveDir 'Z2_190701to190822_YmaxFull_' saveTag '.mat'],'YmaxFullAllDays','-v7.3');
     elseif strcmp(batId,'Zu')
        save([saveDir 'Zu_190704to190820_YmaxFull_' saveTag '.mat'],'YmaxFullAllDays','-v7.3');
     elseif strcmp(batId,'Za')
        save([saveDir 'Za_190524to190530_YmaxFull_' saveTag '.mat'],'YmaxFullAllDays','-v7.3');    
    end
end



