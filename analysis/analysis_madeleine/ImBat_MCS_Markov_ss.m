function [out_mcs_markov_ss] = ImBat_MCS_Markov_2()

num_clust = 2;

% Create vector of all flight cluster IDs 

true_times = flightPaths.AllFlightsMasterTime(flightPaths.flight_starts_idx(:));
[tts,tts_idx]  = sort(true_times(:));
c_s = flightPaths.id(tts_idx);

c_s_2groups = [];
for i=1:size(c_s,1)
    if c_s(i) == 1
        c_s_2groups(i) = 1;
    else
        c_s_2groups(i) = 2;
    end
end
    
c_s_combined = [c_s,c_s_2groups'];

% Get the state transition probabilities
ss = [];
for i=1:num_clust
    i_idx = find(c_s_combined(:,2) == i);
    for j=1:size(i_idx,1)
        if i_idx(j) == size(c_s_combined,1)
            disp('')
        else 
            ss(i,j) = c_s_combined(i_idx(j)+1,2);
        end
    end
end

S = zeros([2,2]); % will contain the counts of clusters transitioned to for each flight cluster
for i=1:2
    b = size(nonzeros(ss(i,:)),1); % b is the number of flights of cluster i
    bnc = ss(i,:); % bnc is the vector of flight clusters that follow each flight of type b
    u = nonzeros(unique(bnc));
    for j=1:2
        if ismember(j,bnc)
            f = sum(bnc==j);
            S(i,j) = f; 
        else
            S(i,j) = 0;
        end
    end
end

% Construct naive probability matrix
Snorm = []; 
for i=1:2
    row = S(i,:);
    Snorm(i,:) = row./sum(row);
end

        
        


end
