close all;
batName = 'Gio';
dateSesh = '200325';
sessionType = 'raw vs deconvolved';
timeZoom = 1400:1900;%1:size(results.S,2);
neuron = [6 10 19 22];

rawVsDecon = figure();
sgtitle([batName ' ' dateSesh ' ' sessionType]);
for i = 1:length(neuron)
subplot(length(neuron),1,i)
plot(results.C_raw(neuron(i),timeZoom),'r')
hold on
plot(results.C(neuron(i),timeZoom),'g')
plot(results.S(neuron(i),timeZoom)*10,'b')
if i == 1
text((timeZoom(end)-timeZoom(1)-100),max(results.C_raw(2,timeZoom))-1,['time ' num2str(timeZoom(1)) ':' num2str(timeZoom(end))]);
end
title(['Neuron ' num2str(neuron(i))]);
end


% subplot(3,1,2)
% plot(results.P.kernel_pars(:,1))
% title('Rising kernel')
% subplot(3,1,3)
% plot(results.P.kernel_pars(:,2))
% title('decaying kernel')

saveas(rawVsDecon,['/Users/periscope/Desktop/rawVsDecon/' batName '_' dateSesh '_' sessionType '.tif']);