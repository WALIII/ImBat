function Markov_Output_34 = ImBat_MCS_Markov_34(flightPaths)


% Function that takes the preprocessed days of only feeder 3 or only
% feeder4 and does seperate markov models 

% First, manually load files spec:
% https://docs.google.com/document/d/1_jYwkI2XeG-4RJJal47YoMkvuga8QrSRvTzrURwCClE/edit

[CombinedROI,ROI_Data] = ImBat_GroupCalcium(ROI_Data,cell_registered_struct,aligned_data_struct);
ROI_Data = ImBat_RepairFlightData(ROI_Data);

% Only have day 3 flights 
f3_days = [1,3,5,7,8,10,12,13,15];
f4_days = [2,4,6,9,11,14];
f3_days_exclude_lows_days = [1,5,7,8,10,12,13,15];
all_days = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
day_select = 3;

if day_select == 3
    ROI_Data_34 = ROI_Data(f3_days);
elseif day_select == 5
    ROI_Data_34 = ROI_Data(f3_days_exclude_lows_days);
elseif day_select == 0
    ROI_Data_34 = ROI_Data;
else
    ROI_Data_34 = ROI_Data(f4_days);
end

flightPaths34 = ImBat_GroupFlights(ROI_Data_34,'mtf',master_track_file,'dist',1.2);

num_clust = size(flightPaths34.clusterIndex,2);
% Create vector of all flight cluster IDs 
true_times = flightPaths34.AllFlightsMasterTime(flightPaths34.flight_starts_idx(:));
[tts,tts_idx]  = sort(true_times(:));
c_s_34 = flightPaths34.id(tts_idx);

% Get the state transition probabilities
nxt_clst = []; % each column is the vector of which flight that column index goes to
prv_clst = [];
for i=1:num_clust
    i_idx = find(c_s_34 == i);
    for j=1:size(i_idx,1)
        if i_idx(j) == size(c_s_34,1)
            prv_clst(i,j) = c_s_34(i_idx(j)-1);
        elseif i_idx(j) == 1
            nxt_clst(i,j) = c_s_34(i_idx(j)+1);
        else
            nxt_clst(i,j) = c_s_34(i_idx(j)+1);
            prv_clst(i,j) = c_s_34(i_idx(j)-1);
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
% FPnorm is the transition probability from y to x

% Construct naive probability matrix
Fnorm34 = []; FPnorm34 = [];
for i=1:num_clust
    row = F(i,:);
    rowP = FP(i,:);
    Fnorm34(i,:) = row./sum(row);
    FPnorm34(i,:) = rowP./sum(rowP);
end

Fnorm34_first5 = Fnorm34(1:5,1:5);
Fnorm34_noone = Fnorm34(2:10,2:10);
Fnorm34_first10 = Fnorm34(1:10,1:10);
FPnorm34_first10 = FPnorm34(1:10,1:10);

figure(); h = heatmap(Fnorm34, 'Colormap', copper); h.title(strcat("All Flight Cluster Transition Probability Matrix; Day ",num2str(day_select)));
figure(); h = heatmap(Fnorm34_first10, 'Colormap', copper);  h.title(strcat("First 10 Flight Cluster Transition Probability Matrix; Day ",num2str(day_select)));

% DTMC of first 10 clusters
mc = dtmc(Fnorm34_first10); 
figure(); 
h = graphplot(mc);
c = h.EdgeColor;
h.EdgeColor = 'k';
mcp = nonzeros(reshape(mc.P,[size(mc.P,1)*size(mc.P,1),1])).*5;
h.LineWidth = mcp;
h.MarkerSize = [size(flightPaths34.id(flightPaths34.id==1),1)/10,
                size(flightPaths34.id(flightPaths34.id==2),1)/10,
                size(flightPaths34.id(flightPaths34.id==3),1)/10,
                size(flightPaths34.id(flightPaths34.id==4),1)/10,
                size(flightPaths34.id(flightPaths34.id==5),1)/10,
                size(flightPaths34.id(flightPaths34.id==6),1)/10,
                size(flightPaths34.id(flightPaths34.id==7),1)/10,
                size(flightPaths34.id(flightPaths34.id==8),1)/10,
                size(flightPaths34.id(flightPaths34.id==9),1)/10,
                size(flightPaths34.id(flightPaths34.id==10),1)/10];
            
            
% DTMC of first clusters that are large
mc = dtmc(Fnorm34_first5); 
figure(); 
h = graphplot(mc);
c = h.EdgeColor;
h.EdgeColor = 'k';
mcp = nonzeros(reshape(mc.P,[size(mc.P,1)*size(mc.P,1),1])).*5;
h.LineWidth = mcp;
h.MarkerSize = [size(flightPaths34.id(flightPaths34.id==1),1)/10,
                size(flightPaths34.id(flightPaths34.id==2),1)/10,
                size(flightPaths34.id(flightPaths34.id==3),1)/10,
                size(flightPaths34.id(flightPaths34.id==4),1)/10,
                size(flightPaths34.id(flightPaths34.id==5),1)/10];
            

% Look at transitions between stereotyped and non-stereotyped
clear c_s_stereo_34;
c_s_stereo_34 = [];
for i=1:size(c_s_34,1)
    if c_s_34(i) == 1
        c_s_stereo_34(i) = 2;
    else%if c_s_34(i) <= 10
        c_s_stereo_34(i) = 1;
    end
end

% Calculate transition proability between 1 (stereotyped and 2
% (understereotyped)
clear c_s_T_34;
c_s_T_34 = zeros(2,2); 
for i=1:size(c_s_34,1)
    if i == size(c_s_34,1)
        disp('')
        break
    end
    if c_s_stereo_34(i) == 1 & c_s_stereo_34(i+1) == 1
        c_s_T_34(1,1) = c_s_T_34(1,1)+1;
    elseif c_s_stereo_34(i) == 1 & c_s_stereo_34(i+1) == 2
        c_s_T_34(1,2) = c_s_T_34(1,2)+1;
    elseif c_s_stereo_34(i) == 2 & c_s_stereo_34(i+1) == 1
        c_s_T_34(2,1) = c_s_T_34(2,1)+1;
    elseif c_s_stereo_34(i) == 2 & c_s_stereo_34(i+1) == 2
        c_s_T_34(2,2) = c_s_T_34(2,2)+1;
    end
end

c_s_Tnorm_34_first10(1,:) = c_s_T_34(1,:)./sum( c_s_T_34(1,:));
c_s_Tnorm_34_first10(2,:) = c_s_T_34(2,:)./sum( c_s_T_34(2,:));

% Netowkr of stereoyped to non
mc = dtmc(c_s_Tnorm_34_first10); 
figure(); 
h = graphplot(mc);
c = h.EdgeColor;
h.EdgeColor = 'k';
mcp = nonzeros(reshape(mc.P,[size(mc.P,1)*size(mc.P,1),1])).*8;
h.LineWidth = mcp;
h.MarkerSize = [(c_s_T_34(1,1)+c_s_T_34(1,2))/10,
                (c_s_T_34(2,1)+c_s_T_34(2,2))/10];
            
disp(c_s_T_34(1,1)+c_s_T_34(1,2));
disp(c_s_T_34(2,1)+c_s_T_34(2,2));

OG_Fnorm = Fnorm34;

% Do shuffle of the sequences (x1000) to see if the probability is
% significant compare to the transitions expected by change (physically
% possible transitions)

% Find transitions that never occur
illegal_sequences_all = [];
illegal_sequences_day34 = [];
for i=1:size(Fnorm34,1)
    zero_idxs = find(Fnorm34(i,:) == 0);
    for j=1:size(zero_idxs,2)
        illegal_sequences_day34 = [illegal_sequences_day34;[i,zero_idxs(j)]];
    end
end

% Find legal transitions instead.. put in dictionary
% Weight each vector with the probability of transition

legal_sequences_day34 = cell(1,size(Fnorm34,1));
weighted_legal_sequences_day34 = cell(1,size(Fnorm34,1));
for i=1:size(Fnorm34,1)
    nonzero_idxs = find(Fnorm34(i,:) ~= 0);
    legal_sequences_day34{i} = nonzero_idxs;
    %for j=1:size(nonzero_idxs,2)
    %    legal_sequences_day34{i} = [legal_sequences_day34{i},nonzero_idxs(j)];
    %end
    w = round(Fnorm34(i,nonzero_idxs)*100);
    for k=1:size(legal_sequences_day34{i},2)
       tt = repmat(legal_sequences_day34{i}(k),1,w(k));
       weighted_legal_sequences_day34{i} = [weighted_legal_sequences_day34{i},tt];
    end
end

% Shuffle the weighted_legal_sequences_day34 arrays
for i=1:size(Fnorm34,1)
    WLS_day34{i} = weighted_legal_sequences_day34{i}(randperm(length(weighted_legal_sequences_day34{i})));
end

for i=1:size(Fnorm34,1)
    LS_day34{i} = legal_sequences_day34{i}(randperm(length(legal_sequences_day34{i})));
end

% Randomly generate sequences from the "legal" sequences
num_shuffle = 1000;
clear S;
S = zeros(num_shuffle,size(c_s_34,1));
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
    S(i,:) = new_seq;
end

% Calculate the T matrix for each S
for m=1:size(S,1)
    SS = S(m,:)';
    clear nxt_clst;
    clear prv_clst;
    nxt_clst = []; % each column is the vector of which flight that column index goes to
    prv_clst = [];
    for i=1:num_clust
        i_idx = find(SS == i);
        if isempty(i_idx)
            prv_clst(i,j) = 0;
            nxt_clst(i,j) = 0;
        else
            for j=1:size(i_idx,1)
                if i_idx(j) == size(SS,1)
                    prv_clst(i,j) = SS(i_idx(j)-1);
                elseif i_idx(j) == 1
                    nxt_clst(i,j) = SS(i_idx(j)+1);
                else
                    nxt_clst(i,j) = SS(i_idx(j)+1);
                    prv_clst(i,j) = SS(i_idx(j)-1);
                end
            end
        end
    end

    clear F FP;
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
    % FPnorm is the transition probability from y to x

    % Construct naive probability matrix
    clear Fnorm34; clear FPnorm34;
    Fnorm34 = []; FPnorm34 = [];
    for i=1:num_clust
        row = F(i,:);
        rowP = FP(i,:);
        Fnorm34(i,:) = row./sum(row);
        FPnorm34(i,:) = rowP./sum(rowP);
    end
    
    %figure(); h = heatmap(Fnorm34, 'Colormap', copper); h.title(strcat("Simulated Sequence #",num2str(m),"; Day ", num2str(day_select)));
    FnormD(:,:,m) = Fnorm34;
end
      
FnormM = nanmean(FnormD,3);
figure(); h = heatmap(FnormM, 'Colormap', copper); h.title(strcat("Mean Simulated Sequence; Day ", num2str(day_select)));

% Make joint probability T by mulilpying each row by the probability of i
Joint_OG_Fnorm = [];
for i=1:size(OG_Fnorm,1)
    p_i = size(c_s_34(c_s_34==i),1)/size(c_s_34,1);
    Joint_OG_Fnorm(i,:) = OG_Fnorm(i,:)*p_i;
end
Joint_OG_Fnorm_vec = reshape(Joint_OG_Fnorm,1,size(Joint_OG_Fnorm,1)^2);

% Do the same for simulated data
% Get prior on each flight
for i=1:size(FnormM,1)
    for j=1:size(S,1)
        Srow = S(j,:);
        S_count(j) = size(Srow(Srow==i),2);
    end
    S_count_mean(i) = mean(S_count);
    P_I = S_count_mean./size(c_s_34,1);
end

Joint_sim_Fnorm = [];
for i=1:size(FnormM,1)
    p_i = size(c_s_34(c_s_34==i),1)/size(c_s_34,1);
    Joint_sim_Fnorm(i,:) = FnormM(i,:)*P_I(i);
end
Joint_sim_Fnorm_vec = reshape(Joint_sim_Fnorm,1,size(FnormM,1)^2);

% K-S Test to see if the join probabilities are sig diff.
clear p h;
[h,p] = kstest2(Joint_sim_Fnorm_vec,Joint_OG_Fnorm_vec);
[h_row1,p_row1] = kstest2(FnormM(1,:),OG_Fnorm(1,:));

%% Do this same shuffle but for stereotyped and non-stereotyped

% Randomly generate sequences from the "legal" sequences. All start with 1.
num_shuffle = 1000;
clear SS;
SS = zeros(num_shuffle,size(c_s_34,1)-1);
for i=1:num_shuffle
    clear new_seq;
    new_seq = [];
    bin_seq = randi([1 2], 1, size(c_s_34,1)-1);
    SS(i,:) = bin_seq;
end
SS_init = repmat(2,1000,1);
SS_I = [SS_init,SS];

% Calculate transition proability between 1 (stereotyped and 2
% (understereotyped)

for i=1:size(SS_I,1)
    SS_row = SS_I(i,:);
    clear SS_T;
    SS_T = zeros(2,2); 
    for j=1:size(SS_row,2)
        if j == size(SS_row,2)
            disp('')
            break
        end
        if SS_row(j) == 1 & SS_row(j+1) == 1
           SS_T(1,1) = SS_T(1,1)+1;
        elseif SS_row(j) == 1 & SS_row(j+1) == 2
            SS_T(1,2) = SS_T(1,2)+1;
        elseif SS_row(j) == 2 & SS_row(j+1) == 1
            SS_T(2,1) = SS_T(2,1)+1;
        elseif SS_row(j) == 2 & SS_row(j+1) == 2
            SS_T(2,2) = SS_T(2,2)+1;
        end
    end
    SS_T_Stack(:,:,i) = SS_T;
end

SS_T_Stack_M = nanmean(SS_T_Stack,3);
SS_FnormM(1,:) = SS_T_Stack_M(1,:)./sum(SS_T_Stack_M(1,:));
SS_FnormM(2,:) = SS_T_Stack_M(2,:)./sum( SS_T_Stack_M(2,:));

% Take joint probability 
Joint_sim_SS_Fnorm = [];
for i=1:size(SS_FnormM,1)
    %p_i = size(c_s_34(c_s_34==i),1)/size(c_s_34,1);
    Joint_sim_SS_Fnorm(i,:) = SS_FnormM(i,:)*0.5
end
Joint_sim_SS_Fnorm_vec = reshape(Joint_sim_SS_Fnorm,1,size(SS_FnormM,1)^2);

% Take joint probability of real data
Joint_SS_Fnorm = [];
for i=1:size(SS_FnormM,1)
    p_i = size(c_s_stereo_34(c_s_stereo_34==i),2)./size(c_s_stereo_34,2);
    Joint_SS_Fnorm(i,:) = c_s_Tnorm_34_first10(i,:).*p_i;
end
Joint_SS_Fnorm_vec = reshape(Joint_SS_Fnorm,1,size(c_s_Tnorm_34_first10,1)^2);


% K-S Test to see if the join probabilities are sig diff.
clear p h;
[hSS,pSS] = kstest2(Joint_sim_SS_Fnorm_vec,Joint_SS_Fnorm_vec);
[h_SS_row1,p_SS_row1] = kstest2(Joint_sim_SS_Fnorm(2,:),Joint_SS_Fnorm(2,:));



end
