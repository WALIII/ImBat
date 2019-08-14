%  using data that is 203x360x1000

% data2 = imresize(data,.5);

gauf = 10; % Gauss smoothing
sm = 5; % smoothing
last_frame = 1000;
% Take the first 100 frames
GGx = double(data2(:,:,1:last_frame));
%G2 = convn(GG, (reshape([1 1 1] / 10, 1, 1, [])), 'same');

% get an estimation of the background
X = imgaussfilt3(double(GGx(:,:,:)),gauf)-1;
% lowpass filter data
G2 = imgaussfilt3(double(GGx(:,:,:)),2);
% G2 = (double(GG).^2-(X.^2))./X;

% subtract background
G2 = G2 -X;

% subtract minumum
Gm = min(G2,[],3)-1;
G2 = G2-Gm;

%smooth the result
G2 = convn(G2, (reshape([1 1 1] / sm, 1, 1, [])), 'same');

% Play a video
figure(); colormap(gray); for i = 1:last_frame; imagesc(G2(:,:,i)*20); pause(0.01); end