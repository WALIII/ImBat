% Rough late-night Spatial information calculation 

% params:
nitter = 1000; % # of itterations for the shuffle

% cocncat the'during' and 'post' activity
GG = cat(1, sp_bnd_act, sp_bnd_act_pst)';

% circ shift data
id=randi(size(GG,2),1,size(GG,1));
GG2=cell2mat(arrayfun(@(x) circshift(GG(x,:),[1 id(x)]),(1:numel(id))','un',0));

% plot the sp matrix, and the shuffle
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
for i = 1: nitter;
    id=randi(size(GG,2),1,size(GG,1));
    GG2=cell2mat(arrayfun(@(x) circshift(GG(x,:),[1 id(x)]),(1:numel(id))','un',0));
    plot(mean(GG2),'b');
    Gx_shuff(i,:) = (mean(GG2)); % shuffled PSTHs
end
plot(mean(GG),'r','LineWidth',4); % plot again so line is on top..

legend('true','shifted');
title('true mean(spikeMatrix) vs circshift');

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
figure();
hold on;
histogram(ShuffSI);
plot([TrueSI TrueSI],[0 30],'LineWidth',10);
title('Shuffled SI vs true SI');
legend('shuffeled','true');
