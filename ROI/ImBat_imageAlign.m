function [imagesAligned] = ImBat_imageAlign(image1,image2,image3,IM1,IM2,IM3)

%IM1 = GG1;
%IM2 = GG2;
%IM3 = GG3;

[optimizer, metric] = imregconfig('multimodal');

optimizer.InitialRadius = 0.005;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 500;
transformType = 'affine';


IM2_aligned = IM2;
IM1_tform = imregtform(image1,image2,transformType,optimizer,metric);
IM1_aligned = imwarp(IM1,IM1_tform,'OutputView',imref2d(size(IM2)));
if ~isempty(IM3)
    %IM1_aligned = imregister(IM1, IM3, transformType, optimizer, metric);
    %IM2_aligned = imregister(IM2, IM3, transformType, optimizer, metric);
    IM3_tform = imregtform(image3,image2,transformType,optimizer,metric);
    IM3_aligned = imwarp(IM3,IM3_tform,'OutputView',imref2d(size(IM2)));    
else
    IM3_tform = [];
    IM3_aligned = IM2;
end

imagesAligned.image1 = image1;
imagesAligned.image2 = image2;
imagesAligned.image3 = image3;
imagesAligned.IM1 = IM1;
imagesAligned.IM2 = IM2;
imagesAligned.IM3 = IM3;
imagesAligned.IM1_aligned = IM1_aligned;
imagesAligned.IM2_aligned = IM2_aligned;
imagesAligned.IM3_aligned = IM3_aligned;
imagesAligned.IM1_tform = IM1_tform;
imagesAligned.IM3_tform = IM3_tform;
imagesAligned.transformType = transformType;
imagesAligned.optimizer = optimizer;