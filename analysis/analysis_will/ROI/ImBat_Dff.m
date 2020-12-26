function [Ymax, Y, maxFig] = ImBat_Dff(Y,varargin);

batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0;
plotting = 0;

filt_rad = 1;
filt_alpha = 1;

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'batname'
            batName=varargin{i+1};
        case 'plotting'
            plotting = varargin{i+1};
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

% filtering
h=fspecial('gaussian',filt_rad,filt_alpha);
Y=imfilter(Y,h,'circular');


% ImBat_Dff
scaling = 8;

%labels for loading and saving data if running independent fromImBat_analyze
if loadFlag == 1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
    load([pwd '/processed/Motion_corrected_Data_DS.mat']);
end

% Make df/f image

% Filter movie
Y = (convn(Y, double(reshape([1 1 1] /10, 1, 1, [])), 'same'));


% 
% hFilt1 = 50; %this is for smoothing the image to determine background for dividing it out to normalize the pixel intensities on all layers
% hFilt2 = 25; %this is for smoothing the actual image before the exponential to show the layers, use 25 for rgb overlay and 50 for the difference heat plots
% h1 = fspecial('disk',hFilt1); %normalize background smoothing
% h2 = fspecial('disk',hFilt2); %smooth presentation image
% scaling = 10; %scale up the ymax to make less pixelated
% psf = fspecial('gaussian', round(2*gSiz), gSig);
% ind_nonzero = (psf(:)>=max(psf(:,1)));
% psf = psf-mean(psf(ind_nonzero));
% psf(~ind_nonzero) = 0;   % only use pixels within the center disk
% Y_med = median(vidData.Y,3);
% Ydff = vidData.Y - Y_med; %subtract med
% Ydff_tFilt = medfilt3(Ydff); %temporal filtering
% Ydff_filt = imfilter(Ydff_tFilt,psf,'symmetric'); 
% YmaxFull = max(Ydff_filt,[],3); %take max
% YmaxFull= imresize(YmaxFull,scaling); %resize to eliminate pixelation
% IM1_raw = YmaxFull;
% IM1doub = imresize(double(IM1_raw),2);
% bground1=imfilter(IM1doub,h1,'replicate');%
% IM1_filt=IM1doub./(bground1+5);
% IM1 = mat2gray(IM1);
% IM1= imfilter(IM1,h2,'replicate');





% Take median of movie
Y_med = median(Y,3);

% Subtract median
Y = double(Y-Y_med);

% take max
Ymax = max(Y,[],3);
Ymax = imresize(Ymax,scaling);

if plotting == 1
    maxFig = figure();
    colormap(gray);
    imagesc(Ymax);
    hold on;
    xticks = get(gca,'xtick');
    yticks = get(gca,'ytick');
    scaling  = 1.1; %1.1um per pixel
    newlabelsX = arrayfun(@(ax) sprintf('%g', scaling * ax), xticks, 'un', 0);
    newlabelsY = arrayfun(@(ay) sprintf('%g', scaling * ay), yticks, 'un', 0);
    set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
    title(['Max Projection dFF: ' batName ' ' dateSesh ' ' sessionType]);
    xlabel('um'); ylabel('um');
    
end






