function ImBat_imageLoadAlign(IM1,IM2,IM3)
% Load in Data:
%folder1 = 'za190514';
%folder2 = 'za190515';
%folder3 = 'za190517';
% 
% for i=1:3
%     cd([pwd '/' folders{1}]) 
% end

% This is your alignment script
[IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(IM1,IM2,IM3);


% This is my script
[a1,b1] = CaBMI_XMASS(IM1,IM2,IM3);
[a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);

figure(); 
subplot(121)
image((a1(:,:,:)))
title('un-aligned')
subplot(122)
image((a2(:,:,:)))
title('aligned')