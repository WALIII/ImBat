close all;
batName = 'Gen';
dateSesh = '200520';
sessionType = 'riseDecayKernal';

ROI.ker1 = results.P.kernel_pars(:,1);
ROI.ker2 = results.P.kernel_pars(:,2);

%% Fit the average AR model
%Sort structure by Ker1
[~,idx]=sort([ROI.ker1]);
%ROI=ROI(idx);
%Params
CNMF_Fs = 6;
g = [mean([ROI.ker1]) mean([ROI.ker2])];
%AR(2) model definition
s = [0 0 1 0 0 0 0 zeros(1,100)];                   %spikes
y = zeros(1,length(s));                             %fluorescence
x = linspace(0, (length(y))/CNMF_Fs, length(y));    %time
for n = 3:length(y)
    y(n) = g(1)*y(n-1)+g(2)*y(n-2)+s(n);
end 
y = y./max(y);  
%Fit with a double exponential
fo = fitoptions('Method','NonlinearLeastSquares','StartPoint', [3, 3, 0.5, 0.05],'Lower', [0.1, 0.5, 0.1, 0.001],'Upper', [10, 10, 2, 2]);
ft = fittype('A*(exp(-(x-x0)/td)-exp(-(x-x0)/tr))','options',fo);
f = fit(x',y',ft);
coeffvals= coeffvalues(f);
riseDecayPlot = figure();
plot(x,y,'LineWidth',3,'color','k');    xlim([-2 Inf]);  ylim([0 Inf]);
hold on;    stem(x,s,'LineWidth',4, 'Marker','none','color','r');   
text(8,0.5,{['Rise: ' num2str(coeffvals(3)) ' s'],['Decay: ' num2str(coeffvals(2)) ' s']});
xlabel('Time(s)');  ylabel('Fluorescence');
title([batName ' ' dateSesh ' ' sessionType]);
figure();
plot(f,x,y);

saveas(riseDecayPlot,['/Users/periscope/Desktop/rawVsDecon/' batName '_' dateSesh '_' sessionType '.tif']);