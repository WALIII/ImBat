function ImBat_RGBroi

% Load in Data:
mat1 = imresize(results10.results.Cn,4);
mat2 = imresize(results12.results.Cn,4);
mat3 = imresize(results13.results.Cn,4);
IM1 = mat2gray(mat1);
IM2 = mat2gray(mat2);
IM3 = mat2gray(mat3);

minRow = min([size(IM1,1),size(IM2,1),size(IM3,1)]);
minCol = min([size(IM1,2),size(IM2,2),size(IM3,2)]);

% Alignment script
[IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));


% Will RGB script
[a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
[a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);

figure(); 
subplot(121)
image((a1(:,:,:)))
title('un-aligned')
subplot(122)
image((a2(:,:,:)))
title('aligned')

figure()
image((a2(:,:,:)))
title('Galileo: 3/10 (R), 3/12 (G), 3/15 (B)')
set(gca,'YDir','normal');

%%
imagesc(imresize(results10.results.Cn,4)); colormap(gray);
imagesc(imresize(results12.results.Cn,4)); colormap(gray);
imagesc(imresize(results15.results.Cn,4)); colormap(gray);

figure();
[imRGB] = cat(3,results10.results.Cn(1:112,1:153),results12.results.Cn(1:112,1:153),results15.results.Cn(1:112,1:153));

C = imfuse(results10.results.Cn(1:112,1:153),results12.results.Cn(1:112,1:153),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);

a = imresize(results10.results.Cn,4);
b = imresize(results12.results.Cn,4);
c = imresize(results15.results.Cn,4);
agray = mat2gray(a);
bgray = mat2gray(b);
cgray = mat2gray(c);
abcRGB = cat(3,agray(1:448,1:612),bgray(1:448,1:612),cgray(1:448,1:612));
imagesc(abcRGB);


