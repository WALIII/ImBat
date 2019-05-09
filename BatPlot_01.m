clear all



mov_listing=dir(fullfile(pwd,'*.mat'));
mov_listing={mov_listing(:).name};

figure();
for i = 1:size(mov_listing,2); %size(mov_listing,2)
    
    
    load(mov_listing{i},'recbuf','fs'),
    WAV{i} = recbuf; 
%     FS_Spectrogram(recbuf,fs);
%      pause(0.1)

end


% Make a template:

% cluster data:

TEMPLATE=fb_spectro_navigate(WAV{1});

[song, song_r, align,Motif_ind]= FS_Premotor_WavSort(WAV,TEMPLATE,fs,30000);

%check data;
figure(); imagesc(song);

% package data
GGG = song_r(:,align:end-75000);

figure(); FS_Spectrogram(GGG(1,:),fs);
% format data..
disp('Formatting data...');
counter = 1;
for i = 1: size(GGG,1)
% WAVcell{1}{counter} = song_r(i,(align-5)/25*48000:(5+align)/25*48000+size(TEMPLATE,1))';
WAVcell{1}{counter} = GGG(i,:)';
counter = counter+ 1;
end


% Warp audio:
[WARPED_TIME, WARPED_audio_true, Index,startT,endT,WARPED_audio] = FS_PreMotor_Warp(WAVcell,TEMPLATE);




counter = 1;
iii = 1;
for ii = 1:size(WARPED_audio{iii},2)
[IMAGE(:,:,ii), T2,F2]= FS_Spectrogram(WARPED_audio{iii}(:,ii),fs);
fOut = WARPED_audio{iii}(:,ii);
audioVect(:,counter) = downsample(tsmovavg(rms(abs(fOut),2),'s',500,1),100);
audioVectT(:,counter) = fOut;
audioVect_UW(:,counter) = downsample(tsmovavg(rms(abs(WARPED_audio_true{1}(:,ii)),2),'s',500,1),100); 
audioVect_UWT(:,counter) = WARPED_audio_true{1}(:,ii);

counter = counter+1;
end




clear counter

% Spectral Density Image
[Gconsensus,F,T] = CY_Get_Consensus(audioVectT,fs);
[Gconsensus2,F,T] = CY_Get_Consensus(audioVect_UWT,fs);


% Generate Figures
figure(); imagesc(flipdim(mean(Gconsensus{1},3),1));
