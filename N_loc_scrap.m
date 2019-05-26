

A2 = full(A);
edit plot_contours
figure(); plot(A2(:,1));
figure(); reshape(A2(:,1),512);
figure(); reshape(A2(:,1),512,512);
figure(); imagesc(reshape(A2(:,1),512,512));
figure();

for i = 1:100;
H(:,:,i) = reshape(A2(:,i),512,512);
end


%%%  Scrap

figure();
for i = 1:400;
H(:,:,i) = reshape(A2(:,i),512,512);
end
M = mean(H,3);
M = mat2gray(M);
figure(); imagesc(M,[0 0.5]);
colormap(gray);
cmap = jet(200);
pcolor = ind2rgb(mat2gray(squeeze(H(:,:,i))),jet(200));
figure(); imagesc(pcolor);
for i = 1: 2;
Q1(:,:,1) = squeeze(H(:,:,1));
Q1(:,:,2) = squeeze(H(:,:,2));
Q1(:,:,3) = squeeze(H(:,:,3));
figure(); Q1 = mat2gray(Q1);
imagesc(Q1,[0 0.0001]);
Q3 =  imfuse(Q1,Q2);
end
for i = 1:400;
Gray = mat2gray(squeeze(H(:,:,i)));
RGB1 = cat(3, Gray, Gray, Gray);  % information stored in intensity
RGB2 = Gray(10);
RGB2(end, end, 3) = 0;  % All information in red channel
GrayIndex = uint8(floor(Gray * 255));
rcmap(:,:,1) = cat(1,(0:0.01:1).^2,zeros(1,101),zeros(1,101))';
rcmap(:,:,2) = cat(1,zeros(1,101),(0:0.01:1).^2,zeros(1,101))';
rcmap(:,:,3) = cat(1,zeros(1,101),zeros(1,101),(0:0.01:1).^2)';
rcmap(:,:,4) = cat(1,(0:0.01:1).^2,(0:0.01:1).^2,zeros(1,101))';
rcmap(:,:,5) = cat(1,zeros(1,101),(0:0.01:1).^2,(0:0.01:1).^2)';
rcmap(:,:,6) = cat(1,(0:0.01:1).^2,zeros(1,101),(0:0.01:1).^2)';
if i ==1
RGB3      = ind2rgb(GrayIndex, squeeze(rcmap(:,:,randi(6))));
else
RGB_t      = ind2rgb(GrayIndex, squeeze(rcmap(:,:,randi(6))));
RGB3 =  plus(RGB_t,RGB3);
end
end
figure();
im = imagesc(RGB3);
im.AlphaData = (M*3);
for i = 1:400;
Gray = mat2gray(squeeze(H(:,:,i)));
RGB1 = cat(3, Gray, Gray, Gray);  % information stored in intensity
RGB2 = Gray(10);
RGB2(end, end, 3) = 0;  % All information in red channel
GrayIndex = uint8(floor(Gray * 255));
rcmap(:,:,1) = cat(1,(0:0.01:1).^2,zeros(1,101),zeros(1,101))';
rcmap(:,:,2) = cat(1,zeros(1,101),(0:0.01:1).^2,zeros(1,101))';
rcmap(:,:,3) = cat(1,zeros(1,101),zeros(1,101),(0:0.01:1).^2)';
rcmap(:,:,4) = cat(1,(0:0.01:1).^2,(0:0.01:1).^2,zeros(1,101))';
rcmap(:,:,5) = cat(1,zeros(1,101),(0:0.01:1).^2,(0:0.01:1).^2)';
rcmap(:,:,6) = cat(1,(0:0.01:1).^2,zeros(1,101),(0:0.01:1).^2)';
if i ==1
RGB3      = ind2rgb(GrayIndex, squeeze(rcmap(:,:,randi(6))));
else
RGB_t      = ind2rgb(GrayIndex, squeeze(rcmap(:,:,randi(6))));
RGB3 =  plus(RGB_t,RGB3);
end
end