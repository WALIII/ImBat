function [audio] = ImBat_ConcatAudio;

% concat audio once in the audio folde


% Get all files in dir.

filePath = cd;

%mic = 1;
fs = 192000;
plotFlag = 0;
firstFileNum = 14;
lastFileNum = length(dir([filePath ,'/','*audio_trial_*']))-1;
if plotFlag == 1
    figure();
end

% Get bat and session name:
sla2find = strfind(filePath,'\'); % get this based on the directory ( windows)
dateSesh = filePath(sla2find(end-1)+1:sla2find(end)-1);%  dateSesh = '200522';
batName = filePath(sla2find(end-2)+1:sla2find(end-1)-1); % batName = 'Gen';


for mic_i = 4%:4
    try
        fileFirst = load([filePath '/_' dateSesh '_audio_trial_' num2str(firstFileNum) '.mat']);
    catch
        try
            fileFirst = load([filePath '/',batName,'_' dateSesh '_audio_trial_' num2str(firstFileNum) '.mat']);
        catch
            fileFirst = load([filePath '/Gen_' dateSesh '_audio_trial_' num2str(firstFileNum) '.mat']);
        end
    end
    ttlConCat = fileFirst.recbuf(:,end);
    audioConCat = fileFirst.recbuf(:,mic_i);
    for file_i = firstFileNum:lastFileNum
        %load the current and next files
        if file_i == firstFileNum
            fileCur = fileFirst;
        else
            fileCur = fileNext;
        end
        try
            fileNext = load([filePath '/_' dateSesh '_audio_trial_' num2str(file_i+1) '.mat']);
        catch
            try
                fileNext = load([filePath '/',batName,'_' dateSesh '_audio_trial_' num2str(firstFileNum) '.mat']);
            catch
                fileNext = load([filePath '/Gen_' dateSesh '_audio_trial_' num2str(firstFileNum) '.mat']);
            end
        end
        fs=fileCur.fs;
        event_ttls_cur = fileCur.recbuf(:,end); %trial data for current file
        audioMicCur = fileCur.recbuf(:,mic_i);
        [R,LTcur,UT,LL,UL] = risetime(event_ttls_cur,fs); %find times of ttl pulses in SECONDS
        
        event_ttls_next = fileNext.recbuf(:,end); %trial data for next file
        audioMicNext = fileNext.recbuf(:,mic_i);
        [R,LTnext,UT,LL,UL] = risetime(event_ttls_next,fs); %find times of ttl pulses in SECONDS
        
        extra_end = (length(event_ttls_cur)- (LTcur(end)*fs));
        extra_start1 = LTnext(1)*fs;
        try
            extra_start2 = LTnext(2)*fs;
        catch
            break %when you reach the last file, if the file contains less than 2 ttl, it can break and will not get confused
        end
        %calculate the amount that needs to be cut off from the next file
        %depending on whether it spans more or less than 3 seconds
        if extra_end + extra_start1 >= 3*fs
            cutOut = round(extra_end+extra_start1-(3*fs));
        elseif extra_end + extra_start1 < 3*fs
            cutOut = round(extra_end+extra_start2-(3*fs));
        end
        %concatenate the ttl and audio streams
        try
        ttlConCat = vertcat(ttlConCat,event_ttls_next(cutOut+1:end));
        audioConCat = vertcat(audioConCat,audioMicNext(cutOut+1:end));
        catch
            disp('End of file');
        end
        
        if plotFlag == 1
            %plot current file ttls
            subplot(4,1,1);
            plot(event_ttls_cur);
            hold on
            for i = 1:length(LTcur)
                plot(LTcur(i)*fs,0,'o')
            end
            title(['File num ' num2str(file_i)]);
            %plot next file ttls
            subplot(4,1,2);
            plot(event_ttls_next);
            hold on
            for i = 1:length(LTnext)
                plot(LTnext(i)*fs,0,'o');
            end
            title(['File num ' num2str(file_i+1)]);
            %plot the next file with the 'proper time' cut out
            subplot(4,1,3);
            plot(event_ttls_next(cutOut+1:end))
            title(['Cut version of file num ' num2str(file_i+1)]);
            %plot the fully concatenated version
            subplot(4,1,4);
            plot(ttlConCat);
            title('All files connected');
            drawnow
            
            pause
            clf
        end
    end
    
    %confirm the ttl are lined up correctly every 3 seconds
    %[R,LTall,UT,LL,UL] = risetime(ttlConCat,fs);
    %figure();
    %plot(LTall,1:length(LTall),'*')
    %title('TTL pulses every 3 seconds');
    %plot the mic trace concatenated based off ttl pulses
    %[R,LTmic1,UT,LL,UL] = risetime(audioConCat);
    %LTmic1 = LTmic1/fs;
    %[b,a] = butter(3,5000/(fs/2),'high');
    %audioFilt = filtfilt(b,a,audioConCat);
    %figure();
    %plot(1:length(audioFilt),audioFilt);
    %title(['Mic ' num2str(mic_i) ' Filered']);
    figure();
    plot(1:length(audioConCat),audioConCat);
    title([batName ' ' dateSesh ' Mic ' num2str(mic_i) ' raw']);
    drawnow
    % hold on;
    % for i = 1:length(LTmic1)
    %     plot(LTmic1(i)*fs,0,'o')
    % end
    %play concatenated sound
    %micObj = audioplayer(audioConCat,fs);
    %play(micObj);
    %pause
end

[audio] = ImBat_formatAudio(audioConCat,ttlConCat,dateSesh,batName,fs);


function [audio] = ImBat_formatAudio(audioConCat,ttlConCat,dateSesh,batName,fs);

% Format audio and align timestamps to flight data

% create time vector, set first ttl to zero
time_vect = 1:length(ttlConCat);
time_vect = time_vect/fs;
ds_factor = 1000;
ttl_ds = downsample(ttlConCat,1000);
[a, b] = findpeaks(abs(ttl_ds),'MinPeakProminence',0.3,'MinPeakDistance',50);

% % test time diff:
% [a2, b2] = findpeaks(ttlConCat,'MinPeakProminence',0.3,'MinPeakDistance',50);

% figure();
% hold on;
% plot(abs(ttl_ds));
% plot(b,a,'*');

offset_t = b(1)/(fs/ds_factor);
time_vect = time_vect-offset_t;

% data
audio.data = audioConCat;
audio.time_vect = time_vect;
% downsampled data
% audio.data = audioConCat_DS;
% audio.time_vect_DS

% metadata
audio.fs = fs;
audio.batName = batName;
audio.session = dateSesh;




