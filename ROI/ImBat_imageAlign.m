function [IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(IM1,IM2,IM3)

%IM1 = GG1;
%IM2 = GG2;
%IM3 = GG3;

[optimizer, metric] = imregconfig('multimodal');

optimizer.InitialRadius = 0.005;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 1000;
transformType = 'similarity';

if ~isempty(IM3)
    IM3_aligned = IM3;
    IM1_aligned = imregister(IM1, IM3, transformType, optimizer, metric);
    IM2_aligned = imregister(IM2, IM3, transformType, optimizer, metric);
else
    IM2_aligned = IM2;
    IM1_aligned = imregister(IM1, IM2, transformType, optimizer, metric);
end