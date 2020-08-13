function [Ymax, Y, maxFig] = ImBat_Dff(Y,varargin);

batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0;
% ImBat_Dff
scaling = 10;
minLimMult = 0; %0 5 min limit multiplier for max projections
maxLimMult = 25; %32 20 max limit multiplier for max projections
%filt_rad = 10;
%filt_alpha = 2;
    gSig = 1; 
    gSiz = 4.5*gSig; 

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
% h=fspecial('gaussian',filt_rad,filt_alpha);
% Y=imfilter(Y,h,'circular');
% Y = (convn(Y, single(reshape([1 1 1] /10, 1, 1, [])), 'same'));
 %make gaussian filter based on normcorre function
    psf = fspecial('gaussian', round(2*gSiz), gSig);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;   % only use pixels within the center disk
    Y = imfilter(Y,psf,'symmetric');

% Take median of movie
Y_min = median(Y,3);

% Subtract median
Y = Y-Y_min;

% take max
Ymax = max(Y,[],3);
Ymax = imresize(Ymax,scaling);
maxFig = figure();
colormap(gray);
%imagesc(Ymax,[0 10]);
imagesc(Ymax);%,[abs(min(min(Ymax)))*minLimMult abs(min(min(Ymax)))*maxLimMult]);% max(max(Ymax))*maxLimMult]);
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







