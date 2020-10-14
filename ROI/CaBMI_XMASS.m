function [RGB1 RGB2] = CaBMI_XMASS(GG1,GG2,GG3,varargin);
days = 3;


% Make sure you take the median for the input matrixes..
HL = [0.20 .80];
T = 1:size(GG1,2);
F = 1:size(GG1,1);
movie = 0;

% Manual inputs
vin=varargin;
for i=1:length(vin)
    if isequal(vin{i},'Direct') % manually inputing a sort order
        ROI=vin{i+1};
        plotROI =1;
    elseif isequal(vin{i},'Indirect')
        ri=vin{i+1};
    elseif isequal(vin{i},'hl')
        HL = vin{i+1};
    end
end


Llim = HL(1);
Hlim = HL(2);
%files   = (files - min(files(:))) / (max(files(:)) - min(files(:)))*255;
% im1(:,:,:,1)=  (GG1 - min(GG1(:))) / (max(GG1(:)) - min(GG1(:)));
% im1(:,:,:,2)=  (GG2 - min(GG2(:))) / (max(GG2(:)) - min(GG2(:)));
% im1(:,:,:,3)=  (GG3 - min(GG3(:))) / (max(GG3(:)) - min(GG3(:)));

im1(:,:,:,1)=  mat2gray(GG1);
im1(:,:,:,2)=  mat2gray(GG2);
if isempty(GG3)
    im1(:,:,:,3)=  mat2gray(GG2);
else
    im1(:,:,:,3)=  mat2gray(GG3);    
end

% Mean subtracted/Normalized
for ii = 1:days
    M = squeeze(im1(:,:,:,ii));
    imT = (M./(median(M,3)+.01));
    im2(:,:,:,ii) = imT- median(imT,3);
end
%im2 = im2 - mean(im2(:,:,:,:),3);
%im2 = mat2gray(im2);

for ii = 1:days % normalize each channel...
    M = squeeze(im2(:,:,:,ii));
    im2(:,:,:,ii) = mat2gray(M);%(M - min(M(:))) / (max(M(:)) - min(M(:)));
end

rsjp = imadjust(im1(:),[Llim ; Hlim]);
rsjp2 = imadjust(im2(:),[Llim ; Hlim]);

    RGB1 = reshape(rsjp,[size(im1,1),size(im1,2),size(im1,3),3]);
    RGB2 = reshape(rsjp2,[size(im2,1),size(im2,2),size(im2,3),3]);

%F = flip(F,1);
%RGB1 = RGB1(size(RGB1,1):-1:1,:,:);

%image(T,F,(RGB1));set(gca,'YDir','normal');


if movie ==1;
    
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    filename = 'testAnimated.gif';
    for i = 1:121;
        % Text assignment
        % Draw plot for y = x.^n
        if i < 10
            image = squeeze(RGB2(:,:,i,:));
            if plotROI ==1;
                for ii = 1:4
                    if ii<3
                        image = insertShape(image,'circle',[mean(ROI.coordinates{ii}(:,1)) mean(ROI.coordinates{ii}(:,2)) 10],'Color','g','LineWidth',5);
                    else
                        image = insertShape(image,'circle',[mean(ROI.coordinates{ii}(:,1)) mean(ROI.coordinates{ii}(:,2)) 10],'Color','r','LineWidth',5);
                    end
                    image = insertText(image,[0 0],'E1','FontSize',18,'BoxColor','black','TextColor','green');
                    image = insertText(image,[0 30],'E2','FontSize',18,'BoxColor','black','TextColor','red');
                    
                end
            end
        elseif i < 10
            image = squeeze(RGB2(:,:,i,:));
            image = insertText(image,[0 0],txt,'FontSize',18,'BoxColor','black','TextColor','white');
        elseif i ==60
            txt = ['HIT!'];
            image = squeeze(RGB2(:,:,i,:));
            image = insertText(image,[0 0],txt,'FontSize',18,'BoxColor','black','TextColor','white');
        else
            txt = ['Time post hit: ',num2str(i-60)];
            image = squeeze(RGB2(:,:,i,:));
            image = insertText(image,[0 0],txt,'FontSize',18,'BoxColor','black','TextColor','white');
        end
        
        
        imshow(squeeze(image));
        axis tight manual % this ensures that getframe() returns a consistent size
        % Capture the plot as an image
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if i == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            
            
        else
            
            imwrite(imind,cm,filename,'gif','DelayTime',0,'WriteMode','append');
        end
    end
    figure();
    hold on
    plot(zscore(squeeze(mean(mean(GG1,1),2))),'r');
    plot(zscore(squeeze(mean(mean(GG2,1),2))),'g');
    if ~isempty(GG3)
        plot(zscore(squeeze(mean(mean(GG3,1),2))),'b');
    end
end

if size(RGB1,3) ==1;
    RGB1 = squeeze(RGB1(:,:,1,:));
    RGB2 = squeeze(RGB2(:,:,1,:));
end