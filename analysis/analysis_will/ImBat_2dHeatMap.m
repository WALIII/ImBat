function [NormRateMat, occupancyMap,C_map]=  ImBat_2dHeatMap(flights,spikes,val,varargin)
% Where are place cells in 2D space?


% User params
filt = 1.2;
ConsistMap = 0; % make a behavioral consistancy map
C_map = [];

nparams=length(varargin);
if mod(nparams,2)>0
    error('Parameters must be specified as parameter/value pairs');
end
% User input
for i=1:2:nparams
    switch lower(varargin{i})
        case 'consist_map' % make a consistancy map
            ConsistMap=varargin{i+1};
        case 'filter' % filer radius
            filt=varargin{i+1};
    end
end


% pre-allocate maps with Nans
NormRateMat = NaN(50,50);
occupancyMap = NaN(50,50);

if isempty(spikes)
    spikes =1;
end

% Reconstruct the cursor
locs = spikes;
pks = ones(length(locs),1);

TData2.cursorA = flights(1,:); %(zscore(roi_ave2.interp_dff(1,:)) + zscore(roi_ave2.interp_dff(2,:)))  - (zscore(roi_ave2.interp_dff(3,:)) + zscore(roi_ave2.interp_dff(4,:)));
TData2.cursorB = flights(2,:); %(zscore(roi_ave2.interp_dff(5,:)) + zscore(roi_ave2.interp_dff(6,:)))  - (zscore(roi_ave2.interp_dff(7,:)) + zscore(roi_ave2.interp_dff(8,:)));

counter = 1;

% index into the 2D point in the direct cells, save the index
Rperm = randi(size(TData2.cursorA,2),1,size(locs,2)); % random timestamps
for i = 1:size(locs,2)
    
    % Bin the location
    for ii = 1:val;
        BinR(1,counter) = TData2.cursorA(:,locs(i));
        BinR(2,counter) = -TData2.cursorB(:,locs(i));
        
        % Randomize...
        rBinR(1,counter) = TData2.cursorA(:,Rperm(i));
        rBinR(2,counter) = -TData2.cursorB(:,Rperm(i));
        
        counter = counter+1;
    end
    
end


% calculate occupancy:
x1 = TData2.cursorA(:,:);
y1 = -TData2.cursorB(:,:);
%h=fspecial('gaussian',filt,filt);
x1 = [ -3000 x1 3000];
y1 = [-3000 y1 3000];
[valuesF, centersF] = hist3([x1(:) y1(:)],[50 50]);
valuesF = valuesF+1;
occ2use = valuesF;
valuesF2 = imgaussfilt(valuesF,filt);

% make a corr matrix:

if ConsistMap ==1
  [C_map] =   ImBat_ConsistMap(x1,y1,centersF);
end
% Make heatplot

%figure(2);
%h=fspecial('gaussian',filt,filt);
try
    x = BinR(1,:);
    y = BinR(2,:);
catch
    disp(' no activity in this bin, skipping..');
    return
end

x = [ -3000 x 3000];
y = [-3000 y 3000];
[values, centers] = hist3([x(:) y(:)],[50 50]);
values(1,1) = 0;
values(50,50) = 0;
%values2=imfilter(values,h,'circular');
values2 = imgaussfilt(values,filt);
%imagesc((values2));
NormRateMat = (values2./valuesF2);
NormRateMat = imgaussfilt(NormRateMat,filt);

occupancyMap = occ2use;
occupancyMap(1,1) = 0;
occupancyMap(50,50) = 0;
occupancyMap = imgaussfilt(occupancyMap,filt);
occupancyMap(occupancyMap>1) = 2;
occupancyMap(occupancyMap<1.1) = NaN;

