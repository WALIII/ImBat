
function [Maps] = ImBat_2DheatMap(LX,LY,LZ);

% Plot 2D Place map field

plotting = 0;


% Make 2D histogram:

% Make heatplot
filt = 10;
h=fspecial('gaussian',filt,filt);

% Data
Dat{1} = LX(:);
Dat{2} = LY(:);
Dat{3} = LZ(:);
% Order
Ord{1} = [1,2];
Ord{2} = [1,3];
Ord{3} = [2,3];
% Title
Tit{1} = 'XY';
Tit{2} = 'XZ';
Tit{3} = 'YZ';

for i = 1:3;
    
data(:,1) = Dat{(Ord{i}(1))};
data(:,2) = Dat{(Ord{i}(2))};
d = 150;                               % width for x and y bins
x = -3000:d:3000;                       % range for x bins
y = x;                               % use same range for y bins, can be different
[values, centers] = hist3(data, 'edges', {x-d/2,y-d/2}); % subtract half bin width to center bins

values2=imfilter(values,h,'circular');
%filter artifacts
filt2 = 5;
h2=fspecial('gaussian',filt2,filt2);
values2=imfilter(values2,h2,'circular');

% indexa = -5:.1:5; imagesc(indexa,indexa,zeros(100,100));
if plotting ==1;
    figure(1);
subplot(1,3,i);
imagesc(centers{2}([1 end]),centers{1}([1 end]),values2);
title([Tit{i}, ' Projection']);
end

Maps{i}.centers = centers;
Maps{i}. values2 = values2;
end

