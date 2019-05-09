function ImBat_Process_Videos()

% Load in all Tiffs, and concatonate them 
mov_listing=dir(fullfile(pwd,'*.tiff'));
mov_listing={mov_listing(:).name};

filenames=mov_listing;


disp('Parsing Audio and Video files');

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for i=1:length(mov_listing)

    [path,file,ext]=fileparts(filenames{i});
if i ==1;
    Yf = read_file(filenames{i});
    Yf = single(Yf);
    [d1,d2,T] = size(Yf);
else
    Yf_temp = read_file(filenames{i});
    Yf_temp = single(Yf_temp);
    Yf = cat(3,Yf,Yf_temp);
    clear Yf_temp
end
end

% remove global and local artifacts
disp('remove artifacts');
Yf =  ImBat_denoise(Yf);


% Save and remove video from ram


% Motion correction


% ROI extraction





