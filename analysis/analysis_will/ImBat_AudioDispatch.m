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

% enter bat dir:
for i = 1:length(dates2use);

cd(bat2use);
f = dir('**');
cd(f(3).name);
% run and save all audio files:

[audio] = ImBat_ConcatAudio;
disp('saving audio data');
save(['audio_data_',f(3).name,'.mat'],'audio','-v7.3');
cd(DIR)
end
% collect audio after saving:
