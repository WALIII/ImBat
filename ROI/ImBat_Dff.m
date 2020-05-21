function [Ymax, Y, maxFig] = ImBat_Dff(Y,varargin);

batName = [];
dateSesh = [];
sessionType = [];
filt_rad = 1;
filt_alpha = 1;
loadFlag = 0;

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'batname'
            batName=varargin{i+1};
        case 'datesesh'
            dateSesh = varargin{i+1};
        case 'sessiontype'
            sessionType = varargin{i+1};
        case 'loadflag'
            loadFlag = varargin{i+1};
        case 'filt_rad'
            filt_rad = varargin{i+1};
            filt_alpha = filt_rad;        
    end
end

% ImBat_Dff
scaling = 4;

%labels for loading and saving data if running independent fromImBat_analyze
if loadFlag == 1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
    processedFolders = dir('processed*');
    processedNewest = sort({processedFolders(end).name});
    processedNewest = char(processedNewest);
    load([pwd processed_Newest '/Motion_corrected_Data_DS.mat']);
end

% Make df/f image

% Filter movie
% filtering
h=fspecial('gaussian',filt_rad,filt_alpha);
Y=imfilter(Y,h,'circular');
Y = (convn(Y, single(reshape([1 1 1] /10, 1, 1, [])), 'same'));

% Take median of movie
Y_med = median(Y,3);

% Subtract median
Y = Y-Y_med;

% take max
Ymax = max(Y,[],3);
Ymax = imresize(Ymax,scaling);
maxFig = figure();
colormap(gray);
imagesc(Ymax);
set(gca,'YDir','normal');
hold on;
%xticks = get(gca,'xtick');
%yticks = get(gca,'ytick');
%scaling_pixel  = 1.1; %1.1um per pixel
%newlabelsX = arrayfun(@(ax) sprintf('%g', scaling_pixel * ax), xticks, 'un', 0);
%newlabelsY = arrayfun(@(ay) sprintf('%g', scaling_pixel * ay), yticks, 'un', 0);
%set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
title(['Max Projection dFF: ' batName ' ' dateSesh ' ' sessionType]);
xlabel('um'); ylabel('um');







