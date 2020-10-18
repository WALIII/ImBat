
function FS_AV_Parse_batch(DIR,varargin)
% FS_AV_Parse

% Parse Data from FreedomScopes
%   Created: 2015/08/02
%   By: WALIII
%   Updated: 2019/05/02 % added better wireless and large file handeling
%   By: WALIII

% FS_AV_parse will do several things:
%
%   1. Seperate Audio and Video, and place them into .mat files in a .mat
%      directory
%   2. Make spectrogram .gif files for extracted audio, for perusing
%      manually.
%   3. Uses the dependancy extractmedia(), an amazing script created by Nathan Perkins

% Run in the Directory of the .mov files. FS_AVparse should be run first:
% FS_AVparse-->FS_TemplateMatch-->FS_FS_Plot_ROI--> BatchDff2----> FS_Image_ROI--> FS_FS_Plot_ROI



mat_dir=[pwd,'/mat'];
error_dir =[pwd,'/error'];
plot_spectrogram =0;
frames_per_file = 100;
aspect_update =0;


if exist(mat_dir,'dir') rmdir(mat_dir,'s'); end
% if exist(gif_dir,'dir') rmdir(gif_dir,'s'); end
%if exist(error_dir,'dir') rmdir(error_dir,'s'); end


vin=varargin;
for i=1:length(vin)
    if isequal(vin{i},'mat_dir') % manually inputing a sort order
        mat_dir=[pwd,vin{i+1}]; % for manual mat_dir naming
    elseif isequal(vin{i},'ri')
        ri=vin{i+1};
    elseif isequal(vin{i},'frames_per_file')
        frames_per_file=vin{i+1};
    end
end


mkdir(mat_dir);
% mkdir(gif_dir);
% mkdir(error_dir);



outlier_flag=0;
if nargin<1 | isempty(DIR), DIR=pwd; end

mov_listing=dir(fullfile(DIR,'*.mov'));
mov_listing={mov_listing(:).name};



filenames=mov_listing;


disp('Parsing Audio and Video files');

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for i=1:length(mov_listing)
    clear file
    [path,file,ext]=fileparts(filenames{i});
    fprintf(1,formatstring,round((i/length(mov_listing))*100));
    FILE = fullfile(DIR,mov_listing{i})
    
    mkdir([mat_dir,'/',file,'_extraction']);
    file = [file,'_extraction/',file];
    
    % Break into many smaller filles:
    file_info = VideoReader(FILE);
    
    vtimes = 0:frames_per_file:file_info.Duration;
    vtimes = cat(2,vtimes,file_info.Duration);
    
    
    for iii = 1:size(vtimes,2)-1;
        
        try
            [a_ts, a, v_ts, v] = extractmedia(FILE,vtimes(iii),vtimes(iii+1));
            
            % Format VIDEO DATA
            [video.height, video.width, video.channels] = size(v{1});
            
            if video.height/video.width == 9/16;
                aspect_update =1;
                % calculate resize factor
                resize_factor = video.height/240;
                resize_factor = 1/resize_factor;
                video.resize_factor = resize_factor;
                resize_factor = 0.5; % this is the true resize
                
            else
                aspect_update =0;
                % calculate resize factor
                resize_factor = video.height/240;
                resize_factor = 1/resize_factor;
                video.resize_factor = resize_factor;
            end
            
            for ii = 1: size(v,1)
                %video.frames(:,:,:,ii) = v{ii}-(noise(ii,:)-min(noise(30:end,:)));
                temp = v{ii}; %-(noise(ii,:)-min(noise(30:end)));
                temp = squeeze(mean(temp,3));
                
                if aspect_update == 1; % update to the proper
                    temp =  imresize(temp,[480 640]);
                end
                video.frames(:,:,ii) = imresize(temp,[240, 320]);
                
                v{ii} = []; % empty the buffer
                temp = [];
            end
            
            
            %%%%%
        catch
            disp(' No audio found in file');
            v1 = VideoReader(FILE);
            k = 1;
            v1.CurrentTime = vtimes(iii);
            while(v1.CurrentTime < vtimes(iii+1))
                %   if v1.CurrentTime <= vtimes(iii+1);
                %break
                % end
                try
                    temp  = readFrame(v1);
                catch
                    disp(' no more readble frames');
                    continue
                end
                if k ==1;
                    
                    % Format VIDEO DATA
                    [video.height, video.width, video.channels] = size(temp);
                    
                    if video.height/video.width == 9/16;
                        aspect_update =1;
                        % calculate resize factor
                        resize_factor = video.height/240;
                        resize_factor = 1/resize_factor;
                        video.resize_factor = resize_factor;
                        resize_factor = 0.5; % this is the true resize
                    else
                        aspect_update =0;
                        % calculate resize factor
                        resize_factor = video.height/240;
                        resize_factor = 1/resize_factor;
                        video.resize_factor = resize_factor;
                    end
                end
                temp = squeeze(mean(temp,3));
                temp = imresize(temp,[240, 320]); % hard code the resize
                video.frames(:,:,k) = temp;
                k = k+1;
            end
            a = 0;
            a_ts = 0;
            v_ts = 0;
            clear k v1;
        end
        %%%%%
        video.times = v_ts;% 0.1703*(day-1);
        video.nrFramesTotal = size(video.frames,3);
        video.FrameRate = 1/mean(diff(v_ts));
        
        
        
        fname = [mat_dir,'/',file,'_',sprintf('%03d', iii),'.tiff'];
        FS_tiff(video.frames,'fname',fname);
        
        %dont save video.frames ( its already a tif...
        video.frames = [];
        
        try
            disp('Performing Gain correction')
            noise = squeeze(mean(mean(squeeze(video.frames(:,[1:30, video.height-30:video.height],3,:)),1),2)); % blue channel
            noise2 = squeeze(mean(mean(squeeze(video.frames([1:30 video.width-30:video.width],:,3,:)),1),2)); % blue channel
            noise = (noise+noise2)/2;
            % sig = squeeze(est(:,:,2,:)); % green channel
            % sig = (squeeze(mean(mean(sig(:,1:80,:),1))));
            video.gain = noise;
            clear noise;
            clear est;
        catch
            video.gain = [];
        end
        
        
        % Format AUDIO DATA
        audio.nrChannels = size(a,2);
        audio.bits = 16;
        audio.nrFrames = length(a);
        audio.data = double(a);
        audio.rate = 48000;
        audio.TotalDurration = audio.nrFrames/48000;
        audio.times = a_ts;
        mic_data = double(a);
        fs = 48000;
        
        
        
        if plot_spectrogram ==1;
            [b,a]=ellip(5,.2,80,[500]/(fs/2),'high');
            plot_data=mic_data./abs(max(mic_data));
            try
                [s,f,t]=fb_pretty_sonogram(filtfilt(b,a,mic_data./abs(max(mic_data))),fs,'low',2.5,'zeropad',0);
                
                
                minpt=1;
                maxpt=min(find(f>=10e3));
                
                imwrite(flipdim(uint8(s(minpt:maxpt,:)),1),hot,fullfile(gif_dir,[file '.gif']),'gif');
            catch
                disp('no audio... skipping spectrogram');
            end
        end
        save(fullfile(mat_dir,[file,'_',sprintf('%03d', iii), '.mat']),'audio','video','-v7.3');
        
        % clear the buffer
        clear video  audio a_ts a v_ts v;
        aspect_update = 0;
    end
end
fprintf(1,'\n');
%%
