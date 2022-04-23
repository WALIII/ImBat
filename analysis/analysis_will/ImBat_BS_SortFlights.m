function [BestScore] = ImBat_BS_SortFlights(flightpathsA,flightpathsB);
% ImBat_BS_SortFlights.m

% make sure that flights that were clustered on different days are the
% same, and generate a confidence score


% To Do: sort corr matrix based on max values


% Get flight data for flightPathA
flight2use = 1:max(flightpathsA.id);
for ii = flight2use
% build flight paths:
for i = 1:size(flightpathsA.flight_starts_idx,2)
    try
        Flights2Use(:,:,i) =  flightpathsA.AllFlights(flightpathsA.flight_starts_idx(i)-120:flightpathsA.flight_starts_idx(i)+(10*120),:);
    catch
        Flights2Use(:,:,i) =  flightpathsA.AllFlights(flightpathsA.flight_starts_idx(i-1)-120:flightpathsA.flight_starts_idx(i-1)+(10*120),:);
    end
end
% create Mean traj
idx2use = find(flightpathsA.id == ii);
MF2U(:,:,ii) = nanmean(Flights2Use(:,:,idx2use),3);
clear Flights2Use;
clear Flights2Use idx2use;
end

% Get flight data for flightPathB
flight2use = 1:max(flightpathsB.id);
for ii = flight2use
% build flight paths:
for i = 1:size(flightpathsB.flight_starts_idx,2)
    try
        Flights2Use(:,:,i) =  flightpathsB.AllFlights(flightpathsB.flight_starts_idx(i)-120:flightpathsB.flight_starts_idx(i)+(10*120),:);
    catch
        Flights2Use(:,:,i) =  flightpathsB.AllFlights(flightpathsB.flight_starts_idx(i-1)-120:flightpathsB.flight_starts_idx(i-1)+(10*120),:);
    end
end
% create Mean traj
idx2use = find(flightpathsB.id == ii);
MF2U_2(:,:,ii) = nanmean(Flights2Use(:,:,idx2use),3);
clear Flights2Use idx2use;
end


% now, run the corr for all dims
for fltsA = 1:size(MF2U,3);
for fltsB = 1:size(MF2U_2,3);
for dims = 1:3;
   score2use(dims) =  corr(MF2U(:,dims,fltsA), MF2U_2(:,dims,fltsB));
end
   score2use = nanmean(score2use);
corrMat(fltsA,fltsB) = score2use;
   clear score2use

end
end

corrMatTemp = corrMat;
for i = 1: size(corrMat,2)
[aa BestScore(i)] = max(corrMatTemp(:,i));
BestScore(i);
if BestScore(i)==1;
    BestScore(i) = tempID;
corrMatTemp(BestScore(i),:) = 0;
corrMatTemp(:,BestScore(i)) = 0;
end

% remove ones with index
Index2use = 1:size(corrMat,2);

unique_early = find(BestScore(2:size(corrMat,1)) ==1);
replaceNums = (length(BestScore)-length(unique_early))+1:length(BestScore)
BestScore(unique_early) = replaceNums;
Index2use(BestScore) = [];
BestScore(find(BestScore(1:end)==1)) = [1 Index2use];
disp('end');% corrMat_itter = corrMat;
% % score
% [aa BestScore] = max(corrMat);
% % assign all extra ones to a unique val
% toFreeze = find(BestScore==1); 
% BestScore(toFreeze) = toFreeze;
% ind2itter = find(hist(BestScore,unique(BestScore))>1);
% ind2lock = find(hist(BestScore,unique(BestScore))==1);
% corrMat_itter(ind2lock,:) = 0;
% corrMat_itter(:,ind2lock) = 0;
% 
% % unique score
% ind2itter = find(hist(BestScore,unique(BestScore))>1)
% ind2itter(ind2itter==1) = []; % just to check
% BestScore_residual = BestScore;
% 
% for i = 1:length(ind2itter);
% checkval = ind2itter(i)
% ind2check = find( BestScore_residual==checkval) 
% [~, ind2use] = max(aa(ind2check))
% % update BestScore, remove this unique relation
% corrMat_itter(ind2check(ind2use),:) = 0;
% corrMat_itter(:,ind2check(ind2use)) = 0;
% ind2check(ind2use) = [];
% [aa BestScore_residual] = max(corrMat_itter);
% end


