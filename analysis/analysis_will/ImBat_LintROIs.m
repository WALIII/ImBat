function score = ImBat_LintROIs(ROI_Data,day2use);

% Get ROI SNR score, and watershed score, then use this to prune ROIs

% PROTOTYPE: 


% WAL3
% 12/14/2020

% 
rmvThresh = 1;
pauseOn = 0;
ii = day2use;% day
A = full(ROI_Data{ii}.ROIs.results.A);
cn = ROI_Data{1, 1}.ROIs.results.Cn;
A = reshape(A,size(cn,1),size(cn,2),size(A,2));


figure();
hold on;
for i =1:size(  ROI_Data{ii}.ROIs.results.C_raw,1);
clf('reset');
subplot(1,5,1:4)
hold on;

% remove slow detrend:
xx = smooth(ROI_Data{ii}.ROIs.results.C_raw(i,:),10);
%gg2 = highpass(gg,0.01,30,'ImpulseResponse','iir','Steepness',0.5);
xx = xx(10000:end-10000);


trend_dat = smooth(xx,10000);
detrended_data = xx-trend_dat;
% detrended_data = xx*1;
% plot(detrended_data,'r')
% xx = smooth(detrended_data,5);
% xx = zscore(xx);
detrended_data = detrended_data-median(detrended_data);
detrended_data = smooth(detrended_data,30);

a1 = prctile(detrended_data,99.99);
a2 = prctile(detrended_data,70);
plot(detrended_data);
detrended_data2(i,:) = detrended_data;
subplot(1,5,5)
imagesc(squeeze(A(:,:,i)));

score(i) = a1-(a2);
title([' ROI: ',num2str(i),' Score: ', num2str(score(i))])

%pause();
end

disp(['Removal threshold set to ',num2str(rmvThresh)]);
% Plot rejected ROIs
ind2rmv = find(score<rmvThresh);
perrmvd = (length(ind2rmv)./size(ROI_Data{ii}.ROIs.results.C_raw,1))*100;

disp([num2str(length(ind2rmv)),' Flaged for Removal, which is ',num2str(perrmvd),' percent of ROIs']);

figure();
hold on;
for i =ind2rmv;
clf('reset');
subplot(1,5,1:4)
hold on;
xx = smooth(ROI_Data{ii}.ROIs.results.C_raw(i,:),10);
xx = zscore(xx);
xx = xx-mode(xx);
%% Spatial criterion
ff = squeeze(A(:,:,i));
[~,aa2] = max(mean(ff));
[~,aa1] = max(mean(ff,2));

if size(ff,2)- aa2 <5 || size(ff,1)- aa1 < 5
    sptit = 'Contains Spatial Anomlies';
    score(i) = score(i)/2; % divide score by two if there are anomolies..
else
   sptit = '';
end
if score(i)<rmvThresh/2
    title('THROW');
else
    title('KEEP');
end

    plot(detrended_data2(i,:));
subplot(1,5,5)
imagesc(squeeze(A(:,:,i)));

%% Spatial criterion
ff = squeeze(A(:,:,i));
[~,aa2] = max(mean(ff));
[~,aa1] = max(mean(ff,2));

if size(ff,2)- aa2 <5 || size(ff,1)- aa1 < 5
    sptit = 'Contains Spatial Anomlies';
    score(i) = score(i)/2; % divide score by two if there are anomolies..
else
   sptit = '';
end




title([' ROI: ',num2str(i),' Score: ', num2str(score(i)), sptit])
if pauseOn ==1;
pause();
end
end
clear ind2rmv perrmvd
ind2rmv = find(score<rmvThresh/2);
perrmvd = (length(ind2rmv)./size(ROI_Data{ii}.ROIs.results.C_raw,1))*100;
disp([num2str(length(ind2rmv)),' ROIs Removed, which is ',num2str(perrmvd),' percent of ROIs']);




% summary stats

figure(); hold on; histogram(score,100); 
line([0.5 0.5],[0 5],'color','r','LineWidth',1)
%% Scrap
% 
% for i =ind2rmv;
% ff = squeeze(A(:,:,i));
% [~,aa2] = max(mean(ff))
% [~,aa1] = max(mean(ff,2))
% 
% if size(ff,2)- aa2 <5 || size(ff,1)- aa1 < 5
%     title('Will be removed');
% else
%     title('Will not be removed');
% end
% % ff = ff(:);
% % ff(ff ==0) = [];
% % i
% % sum(ff)
% % size(ff,1)
% % mean(ff)./size(ff,1)
% 
% pause(); 
% end
