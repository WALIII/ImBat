%function ImBat_test_medFiltAutoregression
Afull = full(results.A); %convert matrix from sparse to full
Amat = reshape(Afull,60,80,[]); %reshape so 4800 is now 60x80 pixels
trace_full = zeros(2,length(Y(1,1,:)));
trace = zeros(2,length(Y(1,1,6:end-5)));
trace_M3 = zeros(2,length(Y(1,1,6:end-5)));
trace_M5 = zeros(2,length(Y(1,1,6:end-5)));
trace_M9 = zeros(2,length(Y(1,1,6:end-5)));
maxIter = 2;

smoothVelocity = zscore(smooth(flightPaths.batSpeed,100));

deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ...% optimize the baseline);
    'max_tau', 100);    % maximum decay time (unit: frame);

for i = 1:2
trace_temp = mean((Amat(:,:,i).*Y),[1,2]); %take mean of each frame within the first cell ID from cnmfe
trace_full(i,:) = reshape(trace_temp,1,[]);
trace(i,:) = trace_full(i,6:end-5);
trace(i,:) = trace(i,:)-min(trace(i,:));
trace_M3(i,:) = medfilt1(trace(i,:),3);
trace_M5(i,:) = medfilt1(trace(i,:),5);
trace_M9(i,:) = medfilt1(trace(i,:),9);
trace_M9_smooth(i,:) = smooth(trace_M9(i,:),2);

deconv_options.type = 'ar1';

[c3_ar1(i,:), s3_ar1(i,:), options] = deconvolveCa(trace_M3(i,:), deconv_options,'maxIter', maxIter);
[c5_ar1(i,:), s5_ar1(i,:), options] = deconvolveCa(trace_M5(i,:), deconv_options,'maxIter', maxIter);
[c9_ar1(i,:), s9_ar1(i,:), options] = deconvolveCa(trace_M9(i,:), deconv_options,'maxIter', maxIter);
[c9_ar1_smooth(i,:), s9_ar1_smooth(i,:), options] = deconvolveCa(trace_M9_smooth(i,:), deconv_options,'maxIter', maxIter);

deconv_options.type = 'ar2';

[c3_ar2(i,:), s3_ar2(i,:), options] = deconvolveCa(trace_M3(i,:), deconv_options,'maxIter', maxIter);
[c5_ar2(i,:), s5_ar2(i,:), options] = deconvolveCa(trace_M5(i,:), deconv_options,'maxIter', maxIter);
[c9_ar2(i,:), s9_ar2(i,:), options] = deconvolveCa(trace_M9(i,:), deconv_options,'maxIter', maxIter);
[c9_ar2_smooth(i,:), s9_ar2_smooth(i,:), options] = deconvolveCa(trace_M9_smooth(i,:), deconv_options,'maxIter', maxIter);

deconv_options.type = 'exp2';

[c3_exp2(i,:), s3_exp2(i,:), options] = deconvolveCa(trace_M3(i,:), deconv_options,'maxIter', maxIter);
[c5_exp2(i,:), s5_exp2(i,:), options] = deconvolveCa(trace_M5(i,:), deconv_options,'maxIter', maxIter);
[c9_exp2(i,:), s9_exp2(i,:), options] = deconvolveCa(trace_M9(i,:), deconv_options,'maxIter', maxIter);
[c9_exp2_smooth(i,:), s9_exp2_smooth(i,:), options] = deconvolveCa(trace_M9_smooth(i,:), deconv_options,'maxIter', maxIter);


end
%%
for p = 1:2
plot_flightsVSar(p) = figure('units','normalized','outerposition',[0 0 0.5 1]);

p1 = subplot(4,1,1); %plot median filter 3
plot(smoothVelocity(21:end-20),'color','r')
hold on
ylabel('velocity (cm/s)');
yt = get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
set(gca,'xticklabel',{[]});
ylim([0 7]);
title('velocity')

p2 = subplot(4,1,2); %plot median filter 3
plot(c3_ar1(p,:))
hold on; 
plot(c3_ar2(1,:))
plot(c3_exp2(1,:))
plot(trace_M3(1,:))
title(['cell #' num2str(p) ': Median filter 3'])
legend('ar1','ar2','exp2','raw')
set(gca,'xticklabel',{[]});

p3 = subplot(4,1,3); %plot median filter 5
plot(c5_ar1(p,:))
hold on;
plot(c5_ar2(1,:))
plot(c5_exp2(1,:))
plot(trace_M5(1,:))
title('Median filter 5')
legend('ar1','ar2','exp2','raw')
set(gca,'xticklabel',{[]});

p4 = subplot(4,1,4); %plot median filter 5
plot(c9_ar1(p,:))
hold on; 
plot(c9_ar2(1,:))
plot(c9_exp2(1,:))
%plot(c9_ar1_smooth(1,:))
%plot(c9_ar2_smooth(1,:))
plot(trace_M9(1,:))
title('Median filter 9')
legend('ar1','ar2','exp2','raw','smooth-ar1','smooth-ar2')
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
xlabel('time (s)');

linkaxes([p2 p3 p4], 'xy')
end





