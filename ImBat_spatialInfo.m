function ImBat_spatialInfo

prob_i = []; lambd_i = []; lambd = [];
prob_i = tt ./ nansum(tt); %occupancy probability for each bin id tt=time per bin
prob_i = prob_i + eps; %add a tiny number to avoid /0
lambd_i = hfr * 1000; %mean firing rate for bin i hfr= spikes/sec (in Hz)
lambd_i = lambd_i + eps; %add a tiny number to avoid /0
%         lambd = mean(lambd_i);%mean firing rate across bins, this is subject to being skewed by really high fr bins
lambd = sum( lambd_i .* prob_i ); %mean fr, this is the same as what michael used in his code: sum( lambd_i . prob_i ) see hippo_analyze_place_fields_Gaussian_smoothing_Michael.m
dist_tuning(ncounter,:) = lambd_i .* prob_i;
all_avgfr = [];
all_avgfr = lambd;
mi(ncounter) = nansum( (prob_i .* (lambd_i / lambd)) .* log2(lambd_i / lambd) ); %same as Nachum paper

%randomize to determine if selectivity is sig
rand_mi = [];
nperms = 100;
for np = 1 : nperms
    shiftfr = circshift(frnofly,round(rand * size(frnofly,2))); %shuffle spike data with respect to position (frnofly = concatenated spike data from flight trials or continuous)
    rand_h = [];
    for nbins = 1 : 30
        rand_h(nbins) = mean(shiftfr(find(bin == nbins)));
    end
    rand_h(isnan(rand_h)) = 0;
    rand_prob_i = []; rand_lambd_i = []; rand_lambd = [];
    rand_prob_i = tt ./ nansum(tt); %occupancy probability for each bin i
    rand_prob_i = rand_prob_i + eps; %add a tiny number to avoid /0
    rand_lambd_i = rand_h * 1000; %mean firing rate for bin i
    %                 rand_lambd_i(isinf(rand_lambd_i)) = nan; %kick out an NANs that sneak in due to the circular shifting (sometimes there is no time on the ends)
    rand_lambd_i = rand_lambd_i + eps; %add a tiny number to avoid /0
    %         lambd = mean(lambd_i);%mean firing rate across bins, this is subject to being skewed by really high fr bins
    rand_lambd = sum( rand_lambd_i . rand_prob_i ); %mean fr, this is the same as what michael used in his code: sum( lambd_i . prob_i ) see hippo_analyze_place_fields_Gaussian_smoothing_Michael.m
    rand_mi(np) = nansum( (rand_prob_i .* (rand_lambd_i / rand_lambd)) .* log2(rand_lambd_i / rand_lambd) );
end
spatial_pval = [];
spatial_pval = size(find(mi(ncounter)<rand_mi),2) / nperms;
allpvals(ncounter) = spatial_pval;