n= 5; %number of neurons to plot
%cn = 3; %cluster number
colors = ['r', 'g', 'b', 'k'];

figure()
%jj = jet(length(rateMap1D.firingRateSmooth{1}(:,1))*10);
for cn = 1:4
%for each neuron
for i = 1:n %length(rateMap1D.firingRateSmooth{1}(:,1))
    subplot(5,(n/5),i)
    rateMap1D.firingRateSmoothNorm{cn}(i,:) = rateMap1D.firingRateSmooth{cn}(i,:)/max(rateMap1D.firingRateSmooth{cn}(i,:));
    plot(1:length(rateMap1D.firingRateSum{cn}(1,:)),rateMap1D.firingRateSmooth{cn}(i,:),colors(cn))%'Color',jj(i*10,:))
    
hold on
title(['n #' num2str(i)])
%pause
%clf
end
end
sgtitle('Neurons 1-50: Zack_190529_fly-1: rgbk')
hold off

figure()
%jj = jet(length(rateMap1D.firingRateSmooth{1}(:,1))*10);
for cn = 1:4
%for each neuron
for i = n+1:2*n %length(rateMap1D.firingRateSmooth{1}(:,1))
    subplot(5,(n/5),i-n)
    rateMap1D.firingRateSmoothNorm{cn}(i,:) = rateMap1D.firingRateSmooth{cn}(i,:)/max(rateMap1D.firingRateSmooth{cn}(i,:));
    plot(1:length(rateMap1D.firingRateSum{cn}(1,:)),rateMap1D.firingRateSmooth{cn}(i,:),colors(cn))%'Color',jj(i*10,:))
    
hold on
title(['n #' num2str(i)])
%pause
%clf
end
end
sgtitle('Neurons 51-100: Zack_190529_fly-1: rgbk')
hold off

figure()
%jj = jet(length(rateMap1D.firingRateSmooth{1}(:,1))*10);
for cn = 1:4
%for each neuron
for i = (2*n)+1:3*n %length(rateMap1D.firingRateSmooth{1}(:,1))
    subplot(5,(n/5),i-(2*n))
    rateMap1D.firingRateSmoothNorm{cn}(i,:) = rateMap1D.firingRateSmooth{cn}(i,:)/max(rateMap1D.firingRateSmooth{cn}(i,:));
    plot(1:length(rateMap1D.firingRateSum{cn}(1,:)),rateMap1D.firingRateSmooth{cn}(i,:),colors(cn))%'Color',jj(i*10,:))
    
hold on
title(['n #' num2str(i)])
%pause
%clf
end
end
sgtitle('Neurons 101-150: Zack_190529_fly-1: rgbk')
hold off