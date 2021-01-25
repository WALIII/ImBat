function ImBat_AudioDispatch

% wrapper function for all audio functions

% run here:
%DIR = Volumes/server_home/users/tobias/flight/audio;

DIR = cd;

% manually entered dates
bat2use = flightPaths.batID;
dates2use = flightPaths.Dates;
% make processed directory:
save_dir = [DIR,'/processed/',bat2use];
mkdir(save_dir);

DIR = [DIR,'/',bat2use];

% enter bat dir:
cd(bat2use);
for i = 1:length(dates2use);

cd(dates2use{i})
f = dir('**');
cd(f(3).name);
% run and save all audio files:

[audio] = ImBat_ConcatAudio;
disp('saving audio data');
save([save_dir,'/audio_data_',f(3).name,'.mat'],'audio','-v7.3');
cd(DIR)
clear audio f
end
% collect audio after saving:
