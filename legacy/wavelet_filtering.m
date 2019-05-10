% Wavelet filtering:


% n = 3;                   % Decomposition Level
% w =  'sym4';              % Near symmetric wavelet
% range = 1:3000;
% WT = wavedec3(GG1(:,:,:),3,w);    % Multilevel 3D wavelet decomposition.
%
%
% A = cell(1,n);
% D = cell(1,n);
% for k = 1:n
%     A{k} = waverec3(WT,'a',k);   % Approximations (low-pass components)
%     D{k} = waverec3(WT,'d',k);   % Details (high-pass components)
% end
%
% % construct pass bands:
% A2 = cat(4,A{2});
% A2 = squeeze(mean(A2,4));
% D2 = cat(4,D{1});
% D2 = squeeze(mean(D2,4));
%
%
% for i = 1:n;
% Ax{i} = mean(A{i},3);
% end
%
% figure();
% for i = range;
%     colormap(gray)
%  subplot(1,4,1)
% imagesc(double(GG1(:,:,i))-Ax{1},[-40 90]);
%
%     subplot(1,4,2)
% imagesc(A{1}(:,:,i)-Ax{1},[-40 90]);
% title('raw video');
%
% subplot(1,4,3)
% imagesc(A{2}(:,:,i)-Ax{1},[-40 90])
% title('wavelet de-noised');
%
%
% pause(0.01);
% end
% 
