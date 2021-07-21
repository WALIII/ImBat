function HMM_output = ImBat_MCS_HMM(flightPaths)

% MCS 3/31/21

% Path to data: /Users/madeleinesnyder/Desktop/Processed/CellReg_files/ROI_Data

% Script for building a Markov Model of the flight behavior.
% First-order transition statistics only in this version
% Based on https://www.mathworks.com/help/stats/hidden-markov-models-hmm.html#f8288

num_clust = size(flightPaths.clusterIndex,2);

% Create vector of all flight cluster IDs 

true_times = flightPaths.AllFlightsMasterTime(flightPaths.flight_starts_idx(:));
[tts,tts_idx]  = sort(true_times(:));
c_s = flightPaths.id(tts_idx);

% State is the flight that came before that flight

nxt_clst = []; % each column is the vector of which flight that column index goes to
prv_clst = [];
for i=1:num_clust
    i_idx = find(c_s(:) == i);
    for j=1:size(i_idx,1)
        if i_idx(j) == size(c_s,1)
            prv_clst(i,j) = c_s(i_idx(j)-1);
        elseif i_idx(j) == 1
            nxt_clst(i,j) = c_s(i_idx(j)+1);
        else
            nxt_clst(i,j) = c_s(i_idx(j)+1);
            prv_clst(i,j) = c_s(i_idx(j)-1);
        end
    end
end

F = zeros([num_clust,num_clust]); % will contain the counts of clusters transitioned to for each flight cluster
FP = zeros([num_clust,num_clust]);
for i=1:num_clust
    b = size(nonzeros(nxt_clst(i,:)),1); % b is the number of flights of cluster i
    bp = size(nonzeros(prv_clst(i,:)),1);
    bnc = nxt_clst(i,:); % bnc is the vector of flight clusters that follow each flight of type b
    bpc = prv_clst(i,:);
    u = nonzeros(unique(bnc));
    up = nonzeros(unique(bpc));
    for j=1:num_clust
        if ismember(j,bnc)
            f = sum(bnc==j);
            F(i,j) = f; 
        else
            F(i,j) = 0;
        end
        if ismember(j,bpc)
            fp = sum(bpc==j);
            FP(i,j) = fp;
        else
            FP(i,j) = 0;
        end
    end
end

% Fnorm is the transition probability of x to y
% FPnorm is the transition probability from y to x - > STATES

% Construct naive probability matrix
Fnorm = []; FPnorm = [];
for i=1:num_clust
    row = F(i,:);
    rowP = FP(:,i);
    Fnorm(i,:) = row./sum(row);
    FPnorm(i,:) = rowP./sum(rowP);
end

% Sequence of outputs
seqs = c_s;

% Sequence of emissions (probability that a possible output (columns) came
% from a given state (row)

% This will be num_clust x num_clust 

% 1-->1     1-->2       1-->3       1-->4

emis = Fnorm;

% Transition probability matrix
trans = FPnorm; 

PSTATES = hmmdecode(seq,trans,emis);

figure(); heatmap(PSTATES(2:end,1:300))

end
