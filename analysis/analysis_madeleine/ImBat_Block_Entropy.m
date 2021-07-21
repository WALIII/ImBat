function [Entropy,Null_Entropy,Null_Max_Entropy] = ImBat_Block_Entropy(Fnorm34,WLS_day34,LS_day34,weighted_legal_sequences_day34,legal_sequences_day34,c_s_34,day_select,num_clust,ngram)

% Calculate the block entropy of a simulated sequnecs using weighted
% transitions

% Shuffle the weighted_legal_sequences_day34 arrays
for i=1:size(Fnorm34,1)
    WLS_day34{i} = weighted_legal_sequences_day34{i}(randperm(length(weighted_legal_sequences_day34{i})));
end

% Randomly generate sequences from the weighted "legal" sequences
num_shuffle = 1000;
clear WS;
WS = zeros(num_shuffle,size(c_s_34,1));
for i=1:num_shuffle
    clear new_seq;
    new_seq = [];
    init = 1;
    for j=1:size(c_s_34,1)
        if j==1
            tvec = WLS_day34{1};
            seq_bloc = tvec(randi(length(tvec)));
            new_seq = [new_seq;seq_bloc];
        else
            tvec = WLS_day34{new_seq(j-1)};
            seq_bloc = tvec(randi(length(tvec)));
            new_seq = [new_seq;seq_bloc];
        end
    end
    WS(i,:) = new_seq;
end

% Calculate the T matrix for each S
for m=1:size(WS,1)
    disp(m)
    SS = WS(m,:)';
    clear nxt_clst;
    clear prv_clst;
    nxt_clst = []; % each column is the vector of which flight that column index goes to
    prv_clst = [];
    for i=1:num_clust
        i_idx = find(SS == i);
        if isempty(i_idx)
            %prv_clst(i,j) = 0;
            nxt_clst(i,j) = 0;
        else
            for j=1:size(i_idx,1)
                if i_idx(j) == size(SS,1)
                    %prv_clst(i,j) = SS(i_idx(j)-1);
                    disp('fix later');
                elseif i_idx(j) == 1
                    nxt_clst(i,j) = SS(i_idx(j)+1);
                else
                    nxt_clst(i,j) = SS(i_idx(j)+1);
                    %prv_clst(i,j) = SS(i_idx(j)-1);
                end
            end
        end
    end

    clear F FP;
    F = zeros([num_clust,num_clust]); % will contain the counts of clusters transitioned to for each flight cluster
    %FP = zeros([num_clust,num_clust]);
    for i=1:num_clust
        b = size(nonzeros(nxt_clst(i,:)),1); % b is the number of flights of cluster i
        %bp = size(nonzeros(prv_clst(i,:)),1);
        bnc = nxt_clst(i,:); % bnc is the vector of flight clusters that follow each flight of type b
        %bpc = prv_clst(i,:);
        u = nonzeros(unique(bnc));
        %up = nonzeros(unique(bpc));
        for j=1:num_clust
            if ismember(j,bnc)
                f = sum(bnc==j);
                F(i,j) = f; 
            else
                F(i,j) = 0;
            end
            %if ismember(j,bpc)
                %fp = sum(bpc==j);
                %FP(i,j) = fp;
            %else
                %FP(i,j) = 0;
            %end
        end
    end

    % Fnorm is the transition probability of x to y
    % FPnorm is the transition probability from y to x

    % Construct naive probability matrix
    clear Fnorm34weight; %clear FPnorm34;
    Fnorm34weight = []; %FPnorm34 = [];
    for i=1:num_clust
        row = F(i,:);
        %rowP = FP(i,:);
        Fnorm34weight(i,:) = row./sum(row);
        %FPnorm34(i,:) = rowP./sum(rowP);
    end
    
    %figure(); h = heatmap(Fnorm34, 'Colormap', copper); h.title(strcat("Simulated Sequence #",num2str(m),"; Day ", num2str(day_select)));
    FnormDweight(:,:,m) = Fnorm34weight;
end
      
FnormMweight = nanmean(FnormDweight,3);
figure(); h = heatmap(FnormMweight, 'Colormap', copper); h.title(strcat("Mean Simulated Sequence; Day ", num2str(day_select)));

%% Take the block entropy of 3-grams

% Create n-gram chunk list of real sequence
seq_words = [];
for i=1:size(c_s_34,1)-ngram+1
    seq_words(i,:) = [c_s_34(i:i+ngram-1)'];
end

% Identify possible permutations 
% For each unique sequence 
permutations = unique(seq_words, 'rows');
perm_count = zeros(1,size(permutations,1));
for i=1:size(permutations,1)
    perm = permutations(i,:);
    for j=1:size(seq_words,1)
        if seq_words(j,:) == perm;
            perm_count(i) = perm_count(i)+1;
        else
            perm_count(i) = perm_count(i);
        end
    end
end

perm_count_pi = perm_count./sum(perm_count);

% Quanitfy the entropy of the real data
for i=1:size(permutations,1)
    E(i) = -(perm_count_pi(i)*log2(perm_count_pi(i)));
end

Entropy = sum(E);

%% Quantify for the null

for k=1:1000
    SEQ = WS(k,:)';
    % Create n-gram chunk list of real sequence
    seq_words = [];
    for i=1:size(SEQ,1)-ngram+1
        seq_words(i,:) = [SEQ(i:i+ngram-1)'];
    end

    % Identify possible permutations 
    % For each unique sequence 
    clear permutations perm_count
    permutations = unique(seq_words, 'rows');
    perm_count = zeros(1,size(permutations,1));
    for i=1:size(permutations,1)
        perm = permutations(i,:);
        for j=1:size(seq_words,1)
            if seq_words(j,:) == perm;
                perm_count(i) = perm_count(i)+1;
            else
                perm_count(i) = perm_count(i);
            end
        end
    end

    perm_count_pi = perm_count./sum(perm_count);

    % Quanitfy the entropy of the real data
    clear E;
    for i=1:size(permutations,1)
        E(i) = -(perm_count_pi(i)*log2(perm_count_pi(i)));
    end

    Null_Entropy(k) = sum(E);
end

%% Max Entropy solution to the problem

% Shuffle the weighted_legal_sequences_day34 arrays
for i=1:size(Fnorm34,1)
    LS_day34{i} = legal_sequences_day34{i}(randperm(length(legal_sequences_day34{i})));
end

% Randomly generate sequences from the "legal" sequences not weighted
% (maxEnt)
num_shuffle = 1000;
clear WS;
WS = zeros(num_shuffle,size(c_s_34,1));
for i=1:num_shuffle
    clear new_seq;
    new_seq = [];
    init = 1;
    for j=1:size(c_s_34,1)
        if j==1
            tvec = LS_day34{1};
            seq_bloc = tvec(randi(length(tvec)));
            new_seq = [new_seq;seq_bloc];
        else
            tvec = LS_day34{new_seq(j-1)};
            seq_bloc = tvec(randi(length(tvec)));
            new_seq = [new_seq;seq_bloc];
        end
    end
    WS(i,:) = new_seq;
end

% Calculate the T matrix for each S
for m=1:size(WS,1)
    SS = WS(m,:)';
    clear nxt_clst;
    clear prv_clst;
    nxt_clst = []; % each column is the vector of which flight that column index goes to
    prv_clst = [];
    for i=1:num_clust
        i_idx = find(SS == i);
        if isempty(i_idx)
            %prv_clst(i,j) = 0;
            nxt_clst(i,j) = 0;
        else
            for j=1:size(i_idx,1)
                if i_idx(j) == size(SS,1)
                    %prv_clst(i,j) = SS(i_idx(j)-1);
                    disp('fix later');
                elseif i_idx(j) == 1
                    nxt_clst(i,j) = SS(i_idx(j)+1);
                else
                    nxt_clst(i,j) = SS(i_idx(j)+1);
                    %prv_clst(i,j) = SS(i_idx(j)-1);
                end
            end
        end
    end

    clear F FP;
    F = zeros([num_clust,num_clust]); % will contain the counts of clusters transitioned to for each flight cluster
    %FP = zeros([num_clust,num_clust]);
    for i=1:num_clust
        b = size(nonzeros(nxt_clst(i,:)),1); % b is the number of flights of cluster i
        %bp = size(nonzeros(prv_clst(i,:)),1);
        bnc = nxt_clst(i,:); % bnc is the vector of flight clusters that follow each flight of type b
        %bpc = prv_clst(i,:);
        u = nonzeros(unique(bnc));
        %up = nonzeros(unique(bpc));
        for j=1:num_clust
            if ismember(j,bnc)
                f = sum(bnc==j);
                F(i,j) = f; 
            else
                F(i,j) = 0;
            end
            %if ismember(j,bpc)
                %fp = sum(bpc==j);
                %FP(i,j) = fp;
            %else
                %FP(i,j) = 0;
            %end
        end
    end

    % Fnorm is the transition probability of x to y
    % FPnorm is the transition probability from y to x

    % Construct naive probability matrix
    clear Fnorm34; clear FPnorm34;
    Fnorm34ME = []; FPnorm34ME = [];
    for i=1:num_clust
        row = F(i,:);
        %rowP = FP(i,:);
        Fnorm34ME(i,:) = row./sum(row);
        %FPnorm34(i,:) = rowP./sum(rowP);
    end
    
    %figure(); h = heatmap(Fnorm34, 'Colormap', copper); h.title(strcat("Simulated Sequence #",num2str(m),"; Day ", num2str(day_select)));
    FnormDME(:,:,m) = Fnorm34ME;
end
      
FnormMME = nanmean(FnormDME,3);
figure(); h = heatmap(FnormMME, 'Colormap', copper); h.title(strcat("Mean Simulated Sequence; Day ", num2str(day_select)));

%% Quantify for the null

for k=1:1000
    SEQ = WS(k,:)';
    % Create n-gram chunk list of real sequence
    seq_words = [];
    for i=1:size(SEQ,1)-ngram+1
        seq_words(i,:) = [SEQ(i:i+ngram-1)'];
    end

    % Identify possible permutations 
    % For each unique sequence 
    clear permutations perm_count
    permutations = unique(seq_words, 'rows');
    perm_count = zeros(1,size(permutations,1));
    for i=1:size(permutations,1)
        perm = permutations(i,:);
        for j=1:size(seq_words,1)
            if seq_words(j,:) == perm;
                perm_count(i) = perm_count(i)+1;
            else
                perm_count(i) = perm_count(i);
            end
        end
    end

    perm_count_pi = perm_count./sum(perm_count);

    % Quanitfy the entropy of the real data
    clear E;
    for i=1:size(permutations,1)
        E(i) = -(perm_count_pi(i)*log2(perm_count_pi(i)));
    end

    Null_Max_Entropy(k) = sum(E);
end

    




