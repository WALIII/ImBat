function [sig, pval, ShuffSI, TrueSI] = fTS_SI(GG)
% Rough late-night Spatial information calculation 


% GG = squeeze(FlightAlignedROI{1}.S(1,1:400,:))';


% params:
nitter = 100; % # of itterations for the shuffle
alphaVal = 0.01;
plotting = 0;
% cocncat the'during' and 'post' activity
%GG = cat(1, sp_bnd_act, sp_bnd_act_pst)';
if size(GG,1)<2
    %disp('Not enough data');
    sig = 0;
    pval = NaN
    ShuffSI = NaN;
    TrueSI = NaN;
    return
end
% check is there is data here:

% circ shift data
id=randi(size(GG,2),1,size(GG,1));
GG2=cell2mat(arrayfun(@(x) circshift(GG(x,:),[1 id(x)]),(1:numel(id))','un',0));

% plot the sp matrix, and the shuffle
if plotting ==1;
figure();
subplot(1,2,1);
imagesc(GG);
title('true spike matrix');
subplot(1,2,2);
imagesc(GG2)
title('circshift spike matrix');


% plot the means of the true, and n itterations of the shuffled data
figure();
hold on;
plot(mean(GG),'r','LineWidth',4);
end

for i = 1: nitter;
    id=randi(size(GG,2),1,size(GG,1));
    GG2=cell2mat(arrayfun(@(x) circshift(GG(x,:),[1 id(x)]),(1:numel(id))','un',0));
    plot(mean(GG2),'b');
    Gx_shuff(i,:) = (mean(GG2)); % shuffled PSTHs
end

if plotting ==1;
plot(mean(GG),'r','LineWidth',4); % plot again so line is on top..

legend('true','shifted');
title('true mean(spikeMatrix) vs circshift');
end
%%  Bastard SI caclulation
Gx_true = mean(GG); %true PSTHs
Gx_shuff_mean = mean(Gx_shuff); %mean shuffled PSTHs ( firing rate proxy

% Calculate total SI for each shuffle
for i = 1: nitter
ShuffSI(i) = nansum((Gx_shuff(i,:)./(Gx_shuff_mean+1e-20)).*log2(Gx_shuff(i,:)./(Gx_shuff_mean+1e-20)));
end

% Calculate true SI
TrueSI = nansum((Gx_true./(Gx_shuff_mean+1e-20)).*log2(Gx_true./(Gx_shuff_mean+1e-20)));


% plot the hist
if plotting ==1;
figure();
hold on;
histogram(ShuffSI);
plot([TrueSI TrueSI],[0 30],'LineWidth',10);
title('Shuffled SI vs true SI');
legend('shuffeled','true');
end

% calculate p
pval = size(find(ShuffSI>TrueSI),2)./size(ShuffSI,2)
if pval<alphaVal
    sig = 1;
else 
    sig = 0;
end

if plotting ==1;
    pause();
end
