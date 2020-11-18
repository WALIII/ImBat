function [figNonrigid,imagesAligned] = ImBat_imageAlign_nonrigid(image1,image2,image3,IM1,IM2,IM3); %  Basic demon registration code. (To easy understand the algorithm)

saveFlag = 1;
% Alpha (noise) constant
alpha=2;

%if you need to compile the c_files
%curDir = pwd;
%cd('C:\Users\tobias\Documents\MATLAB\demon_registration_version_8f\');
% Clean
%clc; clear all; close all;
% Compile the mex files
%compile_c_files

%make saving directory
if saveFlag == 1
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    %saveDir1 = '/Users/periscope/Desktop/analysis/flight/plots/';
    %saveDir1 = 'C:\Users\tobias\Desktop\analysis\plots\';
    if ~exist([saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign'])
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign']);
    else
        disp('You have been working today...');
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign' filesep];
end

% Read two images
%I1=im2double(imread('images/lenag1.png'));  
%I2=im2double(imread('images/lenag2.png')); 

% Set static and moving image to align image 1 to image 2
S=image2; M1=image1;

% Velocity field smoothing kernel
Hsmooth=fspecial('gaussian',[60 60],10);
% The transformation fields
Tx1=zeros(size(M1)); Ty1=zeros(size(M1));
[Sy,Sx] = gradient(S);
for itt=1:200
	    % Difference image between moving and static image
        Idiff=M1-S;
        % Default demon force, (Thirion 1998)
        %Ux = -(Idiff.*Sx)./((Sx.^2+Sy.^2)+Idiff.^2);
        %Uy = -(Idiff.*Sy)./((Sx.^2+Sy.^2)+Idiff.^2);
        % Extended demon force. With forces from the gradients from both
        % moving as static image. (Cachier 1999, He Wang 2005)
        [My,Mx] = gradient(M1);
        Ux = -Idiff.*  ((Sx./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(Mx./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
        Uy = -Idiff.*  ((Sy./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(My./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
 
        % When divided by zero
        Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;
        % Smooth the transformation field
        Uxs=3*imfilter(Ux,Hsmooth);
        Uys=3*imfilter(Uy,Hsmooth);
        % Add the new transformation field to the total transformation field.
        Tx1=Tx1+Uxs;
        Ty1=Ty1+Uys;
        M1=movepixels(IM1,Tx1,Ty1); %make transformation on the original image
end
figNonrigid = figure('units','normalized','outerposition',[0 0 1 0.6]);
sgtitle(['Gal 123 Nonrigid alignment: alpha ' num2str(alpha)]);
ha = tight_subplot(2,7,[.02 .01],[.01 .08],[.01 .01]);
    set(0,'CurrentFigure',figNonrigid);
    axes(ha(1)); imshow(IM1,[]); title('Day 1');
    axes(ha(2)); imshow(image1,[]); title('Day 1 filt');
    axes(ha(3)); imshow(IM2,[]); title('Day 2');
    axes(ha(4)); imshow(image2,[]); title('Day 2 filt');
    axes(ha(5)); imshow(Tx1,[]); title('Transformation Field x');
    axes(ha(6)); imshow(Ty1,[]); title('Transformation Field y');
    axes(ha(7)); imshow(M1,[]); title('Registered Day 1');    
    
% subplot(2,5,1), imshow(IM1,[]); title('Day 1');
% subplot(2,5,2), imshow(IM2,[]); title('Day 2');
% subplot(2,5,3), imshow(M1,[]); title('Registered Day 1');
% subplot(2,5,4), imshow(Tx,[]); title('Transformation Field x');
% subplot(2,5,5), imshow(Ty,[]); title('Transformation Field y');

%% align image 3 to image 2
if ~isempty(image3)
% Set static and moving image
S=image2; M3=image3; %use the filtered images

% Velocity field smoothing kernel
Hsmooth=fspecial('gaussian',[60 60],10);
% The transformation fields
Tx3=zeros(size(M3)); Ty3=zeros(size(M3));
[Sy,Sx] = gradient(S);
for itt=1:200
	    % Difference image between moving and static image
        Idiff=M3-S;
        % Default demon force, (Thirion 1998)
        %Ux = -(Idiff.*Sx)./((Sx.^2+Sy.^2)+Idiff.^2);
        %Uy = -(Idiff.*Sy)./((Sx.^2+Sy.^2)+Idiff.^2);
        % Extended demon force. With forces from the gradients from both
        % moving as static image. (Cachier 1999, He Wang 2005)
        [My,Mx] = gradient(M3);
        Ux = -Idiff.*  ((Sx./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(Mx./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
        Uy = -Idiff.*  ((Sy./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(My./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
 
        % When divided by zero
        Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;
        % Smooth the transformation field
        Uxs=3*imfilter(Ux,Hsmooth);
        Uys=3*imfilter(Uy,Hsmooth);
        % Add the new transformation field to the total transformation field.
        Tx3=Tx3+Uxs;
        Ty3=Ty3+Uys;
        M3=movepixels(IM3,Tx3,Ty3);  %make transformation on the original image
end

    set(0,'CurrentFigure',figNonrigid);
    axes(ha(8)); imshow(IM2,[]); title('Day 2');
    axes(ha(9)); imshow(image2,[]); title('Day 2 filt');
    axes(ha(10)); imshow(IM3,[]); title('Day 3');
    axes(ha(11)); imshow(image3,[]); title('Day 3 filt');
    axes(ha(12)); imshow(Tx3,[]); title('Transformation Field x');
    axes(ha(13)); imshow(Ty3,[]); title('Transformation Field y');
    axes(ha(14)); imshow(M3,[]); title('Registered Day 3'); 

% subplot(2,5,6), imshow(IM2,[]); title('Day 2');
% subplot(2,5,7), imshow(IM3,[]); title('Day 3');
% subplot(2,5,8), imshow(M3,[]); title('Registered day 3');
% subplot(2,5,9), imshow(Tx,[]); title('Transformation x');
% subplot(2,5,10), imshow(Ty,[]); title('Transformation y');
else
    M3 = [];
    Tx3 = [];
    Ty3 = [];
end

if saveFlag == 1
    saveas(figNonrigid,[saveDir filesep 'Gal_200311and20_day1_day2_day3_' datestr(now,'yymmdd-HHMM') 'alignment_nonrigid_alpha' num2str(alpha) '.tif']);
    savefig(figNonrigid,[saveDir filesep 'Gal_200311and20_day1_day2_day3_' datestr(now,'yymmdd-HHMM') 'alignment_nonrigid_alpha' num2str(alpha) '.fig']);
end
%go back to original folder precompiling
%cd(curDir);

imagesAligned.alpha = alpha;
imagesAligned.IM1 = IM1;
imagesAligned.IM2 = IM2;
imagesAligned.IM3 = IM3;
imagesAligned.image1 = image1;
imagesAligned.image2 = image2;
imagesAligned.image3 = image3;
imagesAligned.IM1_aligned = M1;
imagesAligned.IM2_aligned = IM2;
imagesAligned.IM3_aligned = M3;
imagesAligned.Tx1 = Tx1;
imagesAligned.Ty1 = Ty1;
imagesAligned.Tx3 = Tx3;
imagesAligned.Ty3 = Ty3;

%% run rgb overlay, plot, and save
% [a2,b2] = CaBMI_XMASS(IM1,M1,M3);
%     plotTitle = ['nonrigid alpha'  num2str(alpha) ': Gal Full Sesh: Day 1 (r) v 2 (g) v 3 (b)'];
% 
% plotDay123_aligned = figure();
% image((a2(:,:,:)))
% title(plotTitle);
% set(gca,'xticklabel',[],'yticklabel',[]);
% 
% if saveFlag == 1
%         saveas(plotDay123_aligned,[saveDir filesep 'Gal_200311and20_day1_day2_day3_full' datestr(now,'yymmdd-HHMM') 'nonrigid_alpha' num2str(alpha) '.tif']);
%         savefig(plotDay123_aligned,[saveDir filesep 'Gal_200311and20_day1_day2_day3_full' datestr(now,'yymmdd-HHMM') 'nonrigid_alpha' num2str(alpha) '.fig']);
%     end