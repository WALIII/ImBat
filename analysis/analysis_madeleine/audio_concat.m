filePath = [pwd filesep]; %'C:\Users\Tobias\Desktop\analysis\flight\audio_flight\test\song\';
fileParts = fileparts(filePath); Date = fileParts(end-12:end-7);
%mic = 1;

plotFlag = 0;
firstFileNum = 1;
lastFileNum = length(dir([filePath '_' Date '_audio_trial_*']))-1;
if plotFlag == 1
    figure();
end
for mic_i = 1:4
    clear audioConCat; 
    fileFirst = load([filePath '_' Date '_audio_trial_' num2str(firstFileNum)]);
    ttlConCat = fileFirst.recbuf(:,end);
    audioConCat = fileFirst.recbuf(:,mic_i);
    for file_i = firstFileNum:lastFileNum
        %load the current and next files
        if file_i == firstFileNum
            fileCur = fileFirst;
        else
            fileCur = fileNext;
        end
        fileNext = load([filePath '_' Date '_audio_trial_' num2str(file_i+1) '.mat']);
        
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
        
        if size(LTnext) > 8
            continue;
        else
            %concatenate the ttl and audio streams
            ttlConCat = vertcat(ttlConCat,event_ttls_next(cutOut+1:end));
            audioConCat = vertcat(audioConCat,audioMicNext(cutOut+1:end));
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
    [R,LTall,UT,LL,UL] = risetime(ttlConCat,fs);
    figure();
    plot(LTall,1:length(LTall),'*');
    title('TTL pulses every 3 seconds');
    %plot the mic trace concatenated based off ttl pulses
    [R,LTmic1,UT,LL,UL] = risetime(audioConCat);
    LTmic1 = LTmic1/fs;
    [b,a] = butter(3,5000/(fs/2),'high');
    audioFilt = filtfilt(b,a,audioConCat);
%     figure();
%     plot(1:length(audioFilt),audioFilt);
%     title(['Mic ' num2str(mic_i) ' Filered']);
%     figure();
%     plot(1:length(audioConCat),audioConCat);
%     title(['Mic ' num2str(mic_i) ' raw']);
%     drawnow
    % hold on;
    % for i = 1:length(LTmic1)
    %     plot(LTmic1(i)*fs,0,'o')
    % end
    %play concatenated sound
    %micObj = audioplayer(audioConCat,fs);
    %play(micObj);
    %pause
    save(strcat('audioConCat_',num2str(mic_i)),'audioConCat');
end

% Pass the concatenated data to 
% [out, metrics] = ImBat_MCS_alignTimeStamps(audio,video,TS,Markers,mic_data)
% Audio data should be same length as mic data?





