GG = results;
clear IM
counter = 1;
for ii = 1: 100;
    
IMAGE_G = mat2gray(reshape(full(GG.A(:,ii)),60,80));

testIM = imbinarize(IMAGE_G,.0001);
if sum(testIM(:))<200
for i = 1:3    
sf = randi(10,1,1)*0.1;
IM(:,:,counter,i)=IMAGE_G*sf;
end
counter = counter+1;
else
end


end




IM = mat2gray(IM(10:end,:));
IM = squeeze(max(IM,[],3));
% IM = squeeze(max(IM,[],3));
 IM = imresize(IM,2);
%figure();
clear IM
counter = 1;
for ii = 1: 150;
    
IMAGE_G = mat2gray(reshape(full(GG.A(:,ii)),60,80));

testIM = imbinarize(IMAGE_G,.0001);
if sum(testIM(:))<220
for i = 1:3    
sf = randi(10,1,1)*0.3;
IM(:,:,counter,i)=IMAGE_G*sf;
end
counter = counter+1;
else
end


end




IM = mat2gray(IM(10:end,:,:,:));
IM = squeeze(max(IM,[],3));
% IM = squeeze(max(IM,[],3));
 IM = imresize(IM,8);
 IM = imgaussfilt(IM,2);
figure();
imagesc(IM,[0 1]);