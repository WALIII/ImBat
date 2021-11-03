% Script to count all transitions and plot their occurances over the
% session

c_s_34_seg = c_s_34;

five_percent = round(size(c_s_34,1)*(0.05))
co = 11;
length(c_s_34(c_s_34 > co))
%flightPaths34 = flightPaths

%% Extract every 2-gram possible from the data using sliding window
DGM_all = [];
window_size=1;
for i=1:length(c_s_34_seg)-window_size
    two_gram = c_s_34_seg(i:i+window_size);
    DGM_all(i,:) = two_gram;
end

% Remove all 2-grams that bridge sessions
DGM_blacklist = [];
[ss,rr] = sort(flightPaths34.flight_starts_idx);
day_idx = flightPaths34.day(rr);
dur_idx = flightPaths34.dur(rr);
for i=1:length(DGM_all)-2
    if day_idx(i) ~= day_idx(i+1)
        DGM_blacklist = [DGM_blacklist,[i-1,i,i+1]];
    end
end
DGM_pruned_1 = DGM_all;
DGM_pruned_1(DGM_blacklist,:) = [];

% Remove all 2-grams that contain flights greater than the cutoff 
DGM_pruned_2 = [];
for i=1:length(DGM_pruned_1)
    if sum(DGM_pruned_1(i,:) > co) >= 1
        continue
    else
        DGM_pruned_2 = [DGM_pruned_2;DGM_pruned_1(i,:)];
    end
end

% Remove all 3-grams containing a 1
DGM_pruned_3 = [];
for i=1:length(DGM_pruned_2)
    if sum(DGM_pruned_2(i,:)==1) >= 1 
        continue
    else
        DGM_pruned_3 = [DGM_pruned_3;DGM_pruned_2(i,:)];
    end
end

DGM_T_Mat = zeros(length(unique(DGM_pruned_3)),length(unique(DGM_pruned_3)));
for i=1:length(unique(DGM_pruned_3))
    for j=1:length(unique(DGM_pruned_3))
        for m=1:length(DGM_pruned_3)
            if DGM_pruned_3(m,1) == i & DGM_pruned_3(m,2) == j
                DGM_T_Mat(i,j) = DGM_T_Mat(i,j) + 1;
            end
        end
    end
end
figure(); heatmap(DGM_T_Mat)
DGM_T_Mat_Vec = reshape(DGM_T_Mat,length(DGM_T_Mat)*length(DGM_T_Mat),1);

number_of_power_granting_transition_types = length(DGM_T_Mat_Vec(DGM_T_Mat_Vec>20));
power_granting_transition_types = find(DGM_T_Mat_Vec>50);
figure();hold on; bar(DGM_T_Mat_Vec); yline(50);
    

%% Extract every 3-gram possible from the data using sliding window
c_s_34_seg = c_s_34
TGM_all = [];
window_size=2;
for i=1:length(c_s_34_seg)-window_size
    three_gram = c_s_34_seg(i:i+window_size);
    TGM_all(i,:) = three_gram;
end

% Calculate co based on which flight types make up >5% of the data
five_percent = round(size(c_s_34,1)*(0.05));
colormap = distinguishable_colors(co);
length(c_s_34(c_s_34 > co))

% Remove all 3-grams that bridge sessions
TGM_blacklist = [];
[ss,rr] = sort(flightPaths34.flight_starts_idx);
day_idx = flightPaths34.day(rr);
dur_idx = flightPaths34.dur(rr);
for i=1:length(TGM_all)-2
    if day_idx(i) ~= day_idx(i+1)
        TGM_blacklist = [TGM_blacklist,[i-1,i,i+1]];
    end
end
TGM_pruned_1 = TGM_all;
TGM_pruned_1(TGM_blacklist,:) = [];

% Remove all 3-grams that contain flights greater than the cutoff 
TGM_pruned_2 = [];
for i=1:length(TGM_pruned_1)
    if sum(TGM_pruned_1(i,:) > co) >= 1
        continue
    else
        TGM_pruned_2 = [TGM_pruned_2;TGM_pruned_1(i,:)];
    end
end

% Remove all 3-grams containing a 1
TGM_pruned_3 = [];
for i=1:length(TGM_pruned_2)
    if sum(TGM_pruned_2(i,:)==1) >= 1 
        continue;
    else
        TGM_pruned_3 = [TGM_pruned_3;TGM_pruned_2(i,:)];
    end
end

% Calculate 3-count from the TGM matrix data
TGM_T_MATRIX = zeros(co*co,co);
for i=1:co  
    for j=1:co 
        for m=1:length(TGM_pruned_3)                          
            for k=1:co                                  
                if TGM_pruned_3(m,3) == j & TGM_pruned_3(m,2)== i & TGM_pruned_3(m,1) == k    
                    TGM_T_MATRIX((i-1)*co + j,k) = TGM_T_MATRIX((i-1)*co + j,k) + 1;                 
                end
            end
        end
    end
end
    
figure(); heatmap(TGM_T_MATRIX);
TGM_T_MATRIX_Vec = reshape(TGM_T_MATRIX,length(TGM_T_MATRIX)*size(TGM_T_MATRIX,2),1);

number_of_power_granting_transition_types_3 = length(TGM_T_MATRIX_Vec(TGM_T_MATRIX_Vec>20));
power_granting_transition_types_3 = find(TGM_T_MATRIX_Vec>50);
figure(); hold on; bar(TGM_T_MATRIX_Vec); yline(50); 
 
%% For a given cluster, split the flight sequeunces surrounding that cluster (
% in this case, TGM) into groups according to the preceeding flight. 
cluster_number = 2; 
dl = [];
for i=1:length(TGM_pruned_2)
    if TGM_pruned_2(i,3) ~= cluster_number 
        dl = [dl;i];
    end
end
indexes_of_2 = []; indexes_of_4 = []; indexes_of_5 = [];
for i=1:length(TGM_pruned_3)
    if TGM_pruned_3(i,1)== 2 & TGM_pruned_3(i,2)==2 & TGM_pruned_3(i,3) == 2
        indexes_of_2 = [indexes_of_2;i];
    elseif TGM_pruned_3(i,1) == 4 & TGM_pruned_3(i,2)==2 & TGM_pruned_3(i,3) == 2
        indexes_of_4 = [indexes_of_4;i];    
    elseif TGM_pruned_3(i,1) == 5 & TGM_pruned_3(i,2)==2 & TGM_pruned_3(i,3) == 2
        indexes_of_5 = [indexes_of_5;i];
    end
end

% Kernal density plot
figure(); hold on;
pdSix = fitdist(indexes_of_2,'Kernel','BandWidth',4);
x = 0:.1:45;
ySix = pdf(pdSix,x);
plot(x,ySix,'k-','LineWidth',2,'Color','b');

pdSix = fitdist(indexes_of_4,'Kernel','BandWidth',4);
x = 0:.1:45;
ySix = pdf(pdSix,x);
plot(x,ySix,'k-','LineWidth',2,'Color','r');

pdSix = fitdist(indexes_of_5,'Kernel','BandWidth',4);
x = 0:.1:45;
ySix = pdf(pdSix,x);
plot(x,ySix,'k-','LineWidth',2,'Color','g');

% Timeline Line Plot
figure(); hold on; title(strcat("Red 2 2 2; Green 4 2 2; Blue 5 2 2"));
for i=1:length(c_s_34)
    if ismember(i,indexes_of_2)
        line([i i],[4 6],'Color','r','LineWidth',4);
    elseif ismember(i,indexes_of_4)
        line([i i],[2 4],'Color','g','LineWidth',4);
    end
end
for i=1:length(c_s_34)
    if ismember(i,indexes_of_5)
        line([i i],[0 2],'Color','b','LineWidth',4);
    end
end

% Kernal plot of evolution of flights and 1's
%% Find and plot the evolution of transitions from having 1's between them to not having 1's between them.
pair_item_1 = 30 
pair_item_2 = 32

arbitrary = 50;
indexes_1 = []; indexes_2 = []; indexes_3 = [];
for i=1:length(c_s_34)
    if c_s_34(i) == pair_item_1
        ctr = 0;
        for j=1:length(c_s_34)-i
            if c_s_34(i+j) ~=1 & c_s_34(i+j) ~=pair_item_1 & c_s_34(i+j) ~=pair_item_2
                ctr=0;
                break
            elseif c_s_34(i+j) == 1
                ctr = ctr+1;
            elseif c_s_34(i+j) == pair_item_1
%                 ctr=0;
%                 indexes_2 = [indexes_2,i+j];
%                 break
                if ctr ~= 0
                    indexes_1 = [indexes_1,i+1:i+j-1];
                    ctr=0;
                end
                indexes_2 = [indexes_2,i+j];
                break
            elseif c_s_34(i+j) == pair_item_2
                if ctr ~= 0
                    indexes_1 = [indexes_1,i+1:i+j-1];
                    ctr=0;
                end
                indexes_3 = [indexes_3,i+j];  
                break
            end
        end
    elseif c_s_34(i) == pair_item_2
        for j=1:arbitrary
            if c_s_34(i+j) ~=1 & c_s_34(i+j) ~=pair_item_1 & c_s_34(i+j) ~=pair_item_2
                ctr=0;
                break
            elseif c_s_34(i+j) == 1
                ctr=ctr+1;
            elseif c_s_34(i+j) == pair_item_1
                if ctr ~= 0
                    indexes_1 = [indexes_1,i+1:i+j-1];
                    ctr=0;
                end
                indexes_2 = [indexes_2,i+j];
                break
            elseif c_s_34(i+j) == pair_item_2
%                 ctr=0;
%                 indexes_3 = [indexes_3,i+j];  
%                 break
                if ctr ~= 0
                    indexes_1 = [indexes_1,i+1:i+j-1];
                    ctr=0;
                end
                indexes_3 = [indexes_3,i+j];  
                break
            end
        end
    end
end

% Timeline Line Plot
figure(); hold on; title(strcat("Red 1; Green ",num2str(pair_item_1),";"," ","Blue "," ",num2str(pair_item_2))); 
for i=1:length(c_s_34)
    if ismember(i,indexes_1)
        line([i i],[4 6],'Color','r','LineWidth',4);
    elseif ismember(i,indexes_2)
        line([i i],[2 4],'Color','g','LineWidth',4);
    elseif ismember(i,indexes_3)
        line([i i],[0 2],'Color','b','LineWidth',4);
    end
end

