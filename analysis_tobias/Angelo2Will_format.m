function Angelo2Will_format
% Wrapper to convert to the old format...

dateDir = dir('20*');
for i = 1:length(dateDir)
    cd(dateDir(i).name);

%make new folders 
processed_FN = ['processed_',datestr(now,'yyyy_mm_dd__hhMM')];
trackDir = dir('*track.mat');
fileName = extractBefore(trackDir.name,'_track.mat');
extractFolder = [fileName '_extraction'];
mkdir(['extracted' filesep extractFolder filesep processed_FN]);


% copy over tracking file
disp('Copying Track File')
copyfile(trackDir.name, 'extracted');

% copy over alignment file
disp('Copying Alignment File')
copyfile('Alignment.mat', ['extracted' filesep extractFolder filesep processed_FN]);

% convert ROI data
disp('Converting ROI data...');
S = dir(fullfile(pwd,'*.mat'));
N = {S.name};
X = contains(N,'CNMFe');
out = load(N{X}, 'neuron_1');
results = struct(out.neuron_1);
save(['extracted' filesep extractFolder filesep processed_FN filesep 'results'],'results');

% save tif to mat
disp('Loading video...');
Y = read_file('DSampled.tif');
Ysiz = size(Y);
disp('Saving video...');
save(['extracted' filesep extractFolder filesep processed_FN filesep 'Motion_corrected_Data_DS.mat'],'Y','Ysiz')
clear Y Ysiz % clear data for memory resons

cd(dateDir(i).folder);
end
