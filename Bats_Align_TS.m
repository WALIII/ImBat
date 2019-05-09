
function [out] = Bats_Align_TS(audioIn,TS,Markers);

TS = TS(:,2);
audio = audioIn(:,1);
audio = downsample(audio,100);
fs = 48000;

% resample timestamps:
TS = resample(TS,round(48000/120),1);
TS2 = downsample(TS,100);

HH = zscore(smooth(abs(audio)))-min(zscore(smooth(abs(audio))));
HH2 = zscore(smooth(TS2,10));

[Apks,Alocs] = findpeaks(HH,'MinPeakProminence',4,'MinPeakDistance',60);
[Bpks,Blocs] = findpeaks(HH2,'MinPeakProminence',1,'MinPeakDistance',6);

% make better sig
HH1b = zeros(1,length(HH))';
HH1b(Alocs) = 1;

HH2b = zeros(1,length(HH2))';
HH2b(Blocs) = 1;

% offset:
offsetA = Alocs(1);
offsetB = Blocs(1);

%ur = resample(u,3,2);

% Align to first peak
figure(); 
hold on;
plot((1:length(HH))-offsetA,HH,'r');
plot((1:length(HH2))-offsetB,HH2,'b');



out.WarpedMarkers = [];
out.scaleFactor = [];



figure();
hold on;
plot((1:length(HH))-offsetA,HH1b,'r');
%plot((1:length(HH2))-ones(length(Alocs),1),'r*');
plot((1:length(HH2))-offsetB,HH2b,'b');
% plot(Blocs-offsetB,ones(length(Blocs),1),'b*');


