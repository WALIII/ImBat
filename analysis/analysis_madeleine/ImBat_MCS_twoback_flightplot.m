% Script to test if: 
% P(F(t+1), F(t)) == P(F(t+1), F(t), F(t-1))

% Select the portion of the c_s_34 list you want
day_stretch_1 = [12,28]; % non-stationary
day_stretch_2 = [33,49]; % stationary
[ss,rr] = sort(flightPaths34.flight_starts_idx);
day_idx = flightPaths34.day(rr);
start_1 = find(day_idx == day_stretch_2(1));
end_1 = find(day_idx == day_stretch_2(2));
c_s_34_seg = c_s_34(start_1(1):end_1(end));
li=121;


% % Relabel the numbers to plot?
% c_s_34_seg_mapping = {};
% unique_temp = unique(c_s_34_seg);
% for i=1:length(unique_temp)
%     c_s_34_seg_mapping{i} = unique_temp(i);
% end;
% 
% for i=1:length(c_s_34_seg)
%     c_s_34_seg_r(i) = find([c_s_34_seg_mapping{:}] == c_s_34_seg(i));
% end

% Use all data if not doing the rescaling above
c_s_34_seg = c_s_34;
c_s_34_seg_r = c_s_34;
    
% color_scheme = [jet(20); jet(length(unique(c_s_34_seg_r))-10)];
% f = 2;
% L = [];
% for i=2:length(c_s_34_seg_r)-1
%     if c_s_34_seg_r(i) == f
%         L = [L;c_s_34_seg_r(i-1),c_s_34_seg_r(i),c_s_34_seg_r(i+1),i];
%     end
% end

% L_sorted = sortrows(L);
% figure(); hold on; 
% ylim([0,length(L_sorted)]); ylabel('Flight'); xlabel('F(t-1)                                 F(t)                               F(t+1)');
% title(strcat('Flight Type ',num2str(f), ' surrounding flights'));
% for i=1:length(L_sorted)
%     line([0,5], [i,i],'Color',[color_scheme(L_sorted(i,1),:)],'LineWidth',3);
%     line([5,5],[0,length(L_sorted)],'LineWidth',2,'Color','black');
%     line([5,10], [i,i],'Color',[color_scheme(L_sorted(i,2),:)],'LineWidth',3);
%     line([10,10],[0,length(L_sorted)],'LineWidth',2,'Color','black');
%     line([10,15], [i,i],'Color',[color_scheme(L_sorted(i,3),:)],'LineWidth',3);
% end
% li = length(unique(c_s_34_seg_r));
% 
% figure(); hold on;
% for i=1:20
%     line([i,i],[0,5],'LineWidth',2,'Color',[color_scheme(i,:)]);
% end

%%
% Extract 1000 3-gram sequences from the data
% Construct the P(F(t) | F(t-1)) vector
% Construct the P(F(t) | F(t-1),F(t-2)) vector
% KS test these two vectors

% Extract 1000 3-gram sequences from the data
% TGM = [];
% for i = 1:3000
%     three_gram_idx = randsample([2:length(c_s_34_seg_r)-1],1);
%     three_gram = [c_s_34_seg_r(three_gram_idx-1),c_s_34_seg_r(three_gram_idx),c_s_34_seg_r(three_gram_idx+1)];
%     TGM(i,:) = three_gram;
% end

% Extract every 3-gram possible from the data using sliding window
TGM = [];
window_size=2;
for i=1:length(c_s_34_seg_r)-window_size
    three_gram = c_s_34_seg_r(i:i+window_size);
    TGM(i,:) = three_gram;
end

clear B_0_s B_0_s_norm B_0_s_scaled B_1_s B_1_s_norm
% Calculate B_0 from the TGM matrix data
B_0_s = zeros(li,li);
for i=1:li
    for m=1:length(TGM)
        for j=1:li
            if TGM(m,2) == i & TGM(m,3) == j
                B_0_s(i,j) = B_0_s(i,j) + 1;
            end
        end
    end
end
B_0_s_scaled = [];
for i=1:size(B_0_s,1)
    row = B_0_s(i,:);
    row_scaled = row./possible_flight_num;
    B_0_s_scaled(i,:) = row_scaled;
end

B_0_s_norm = [];
for i=1:size(B_0_s,1)
    row = B_0_s_scaled(i,:);
    row_sum = row./sum(row);
    B_0_s_norm(i,:) = row_sum;
end
B_0_s_norm(isnan(B_0_s_norm))=0;

% Construct OneBack vector (length 1000) from B_0_s matrix and TGM.
OneBack_s = [];
for i=1:5000
    possible_flight_num = length(nonzeros(B_0_s_norm(TGM(i,2),:)));
    OneBack_s(i) = B_0_s_norm(TGM(i,2),TGM(i,3));%/possible_flight_num;
end

% Calculate B_1 from the TGM matrix data
B_1_s = zeros(li*li,li);
for i=1:li   
    for j=1:li 
        for m=1:length(TGM)                          
            for k=1:li                                  
                if TGM(m,3) == j & TGM(m,2)== i & TGM(m,1) == k    
                    B_1_s((i-1)*li + j,k) = B_1_s((i-1)*li + j,k) + 1;                 
                end
            end
        end
    end
end
B_1_s_norm = zeros(size(B_1_s,1),size(B_1_s,2));
for i=1:size(B_1_s,1)
    row = B_1_s(i,:);
    B_1_s_norm(i,:) = row./sum(row);
end
B_1_s_norm(isnan(B_1_s_norm))=0;

% Construct TwoBack vector (length 1000) from B_1_s matrix and TGM.
TwoBack_s = [];
for i=1:5000
    P_B_1_idx = ((TGM(i,2)-1)*li) + TGM(i,1);
    TwoBack_s(i) = B_1_s_norm(P_B_1_idx,TGM(i,3));
end

% Construct B_0 and B_1 from the data NOT TGM and compare for sanity
% B_0 = P(F(t+1) | F(t)) COUNTS
% Row index is (F(t)), Column index is F(t+1)
B_0 = zeros(li,li);
for i=1:li
    for m=2:length(c_s_34_seg_r)
        for j=1:li
            if c_s_34_seg_r(m-1) == i & c_s_34_seg_r(m) == j
                B_0(i,j) = B_0(i,j) + 1;
            end
        end
    end
end
B_0_norm = [];
for i=1:size(B_0,1)
    row = B_0(i,:)
    row_sum = row./sum(row)
    B_0_norm(i,:) = row_sum;
end
B_0_norm(isnan(B_0_norm))=0;

OneBack = [];
for i=1:5000
    possible_flight_num = length(nonzeros(B_0_norm(TGM(i,2),:)));
    OneBack(i) = B_0_norm(TGM(i,2),TGM(i,3));%/possible_flight_num;
end

B_1 = zeros(li*li,li);
for i=1:li   
    for j=1:li 
        for m=3:length(c_s_34_seg_r)                          
            for k=1:li                                  
                if c_s_34_seg_r(m-2) == j & c_s_34_seg_r(m-1)== i & c_s_34_seg_r(m) == k    
                    B_1((i-1)*li + j,k) = B_1((i-1)*li + j,k) + 1;                 
                end
            end
        end
    end
end
B_1_norm = zeros(size(B_1,1),size(B_1,2));
for i=1:size(B_1,1)
    row = B_1(i,:);
    B_1_norm(i,:) = row./sum(row);
end
B_1_norm(isnan(B_1_norm))=0;

TwoBack = [];
for i=1:5000
    P_B_1_idx = ((TGM(i,2)-1)*li) + TGM(i,1);
    TwoBack(i) = B_1_norm(P_B_1_idx,TGM(i,3));
end

%% Perform rounds of a KS test to see if Oneback_s sample vector and TwoBack_s sample vector are from the same
% distribution
clear h_s p_s HS_2 PS_2 HS_1 PS_1
OneBack_Matrix = [OneBack_s',OneBack']; TwoBack_Matrix = [TwoBack_s',TwoBack'];
for i=1:1000
    indexes1 = randsample([1:size(OneBack,2)],30);
    [h_s_1,p_s_1] = kstest2(OneBack_Matrix(indexes1,1),TwoBack_Matrix(indexes1,1),'Alpha',0.01/1000);
    HS_1(i) = h_s_1;
    PS_1(i) = p_s_1;
    [h_s_2,p_s_2] = kstest2(OneBack_Matrix(indexes1,2),TwoBack_Matrix(indexes1,2),'Alpha',0.01/1000);
    HS_2(i) = h_s_2;
    PS_2(i) = p_s_2;
end
disp(strcat(num2str(sum(HS_1))," "," of 1000 kstests of OneBack_s vs. TwoBack_s were significantly different"));
disp(strcat(num2str(sum(HS_2))," "," of 1000 kstests of OneBack vs. TwoBack were significantly different"));

% Looks like the only sig diff things are the two difference Twoback vectors... 
OBs = round(OneBack_s,4); OB = round(OneBack,4);
figure(); hold on; plot(cumsum(OBs)); plot(cumsum(OB)); title("OneBack v.s. OneBack_s Cumsums");
TBs = round(TwoBack_s,4); TB = round(TwoBack,4);
figure(); hold on; plot(cumsum(TBs)); plot(cumsum(TB)); title("TwoBack v.s. TwoBack_s Cumsums")
figure(); hold on; plot(cumsum(TBs)); plot(cumsum(OBs)); title("OneBack_s v.s. TwoBack_s Cumsums")
figure(); hold on; plot(cumsum(TB)); plot(cumsum(OB)); title("OneBack v.s. TwoBack Cumsums")


%% Perform rounds of a KS test to see if OneBack_S and Oneback are the same 
clear h_s p_s HS_2 PS_2 HS_1 PS_1
OneBack_Matrix = [OneBack_s',OneBack']; TwoBack_Matrix = [TwoBack_s',TwoBack'];
for i=1:1000
    indexes1 = randsample([1:size(OneBack,2)],30);
    indexes2 = randsample([1:size(OneBack,2)],30);
    [h_s_1,p_s_1] = kstest2(OneBack_Matrix(indexes1,1),OneBack_Matrix(indexes1,2),'Alpha',0.01/1000);
    HS_1(i) = h_s_1;
    PS_1(i) = p_s_1;
    [h_s_2,p_s_2] = kstest2(TwoBack_Matrix(indexes2,1),TwoBack_Matrix(indexes2,2),'Alpha',0.01/1000);
    HS_2(i) = h_s_2;
    PS_2(i) = p_s_2;
end
disp(strcat(num2str(sum(HS_1))," "," of 1000 kstests of OneBack_s vs. OneBack were significantly different"));
disp(strcat(num2str(sum(HS_2))," "," of 1000 kstests of TwoBack vs. TwoBack_s were significantly different"));

%% See how different the two T matricies are
B_0_KL = [];
for i=1:size(B_0_norm,1)
    row = B_0_norm(i,:);
    row_s = B_0_s_norm(i,:);
    B_0_KL(i) = getKullbackLeibler(row,row_s)
end
B_1_KL = [];
for i=1:size(B_1_s_norm,1)
    row = B_1_norm(i,:);
    row_s = B_1_s_norm(i,:);
    try
        temp = getKullbackLeibler(row,row_s);
        if isinf(temp)
            temp=1;
        end
        B_1_KL(i) = temp;
    catch
        disp("One row is zero, leave out");
        continue;
    end
end
 
mean(B_0_KL)
mean(B_1_KL)

%% Do a null to see if there is significance
% Sort by middle column
% Shuffle the left-most vector in the TGM matrix
[TGM_sorted,TGM_ss] = sort(TGM(:,2));
temp = TGM(:,1);
TGM_sorted_1 = temp(TGM_ss);
temp = TGM(:,3);
TGM_sorted_2 = temp(TGM_ss);
TGM_sorted_full = [TGM_sorted_1,TGM_sorted,TGM_sorted_2];
temp_unique = unique(TGM(:,2));
clear subsegment_matrix;
for i=1:length(temp_unique)
    subsegs = [];
    for j=1:length(TGM_sorted_full)
        if TGM_sorted_full(j,2) == temp_unique(i)
            subsegs = [subsegs; TGM_sorted_full(j,:)];
        end
    end
    subsegment_matrix{i} = subsegs;
end

TGM_sorted_full_shuff = [];
for i=1:length(temp_unique)    
    temp = subsegment_matrix{i};
    temp_1 = temp(:,1);
    temp_1_shuff = temp_1(randperm(length(temp_1)));
    temp_shuff = [temp_1_shuff,temp(:,2),temp(:,3)];
    TGM_sorted_full_shuff = [TGM_sorted_full_shuff;temp_shuff];
end

%% Now do the test for the shuffled data matrix
f=5;  
TGM_test = []; TGM_test_shuff = [];
for i=1:length(TGM_sorted_full)
    if TGM_sorted_full(i,2) == f
        TGM_test = [TGM_test;TGM_sorted_full(i,:)];
    end
    if TGM_sorted_full_shuff(i,2) == f
        TGM_test_shuff = [TGM_test_shuff;TGM_sorted_full_shuff(i,:)];
    end
end

clear OneBack_f OneBack_f_shuff TwoBack_f TwoBack_f_shuff 
% Construct OneBack and TwoBacks using B_0_norm and B_1_norm
OneBack_f = [];
for i=1:length(TGM_test)
    possible_flight_num = length(nonzeros(B_0_norm(TGM_test(i,2),:)));
    OneBack_f(i) = B_0_norm(TGM_test(i,2),TGM_test(i,3));%/possible_flight_num;
end

OneBack_f_shuff = [];
for i=1:length(TGM_test_shuff)
    possible_flight_num = length(nonzeros(B_0_norm(TGM_test_shuff(i,2),:)));
    OneBack_f_shuff(i) = B_0_norm(TGM_test_shuff(i,2),TGM_test_shuff(i,3));%/possible_flight_num;
end

OneBack_f_s = [];
for i=1:length(TGM_test_shuff)
    possible_flight_num = length(nonzeros(B_0_s_norm(TGM_test(i,2),:)));
    OneBack_f_s(i) = B_0_s_norm(TGM_test(i,2),TGM_test(i,3));%/possible_flight_num;
end

OneBack_f_s_shuff = [];
for i=1:length(TGM_test_shuff)
    possible_flight_num = length(nonzeros(B_0_s_norm(TGM_test_shuff(i,2),:)));
    OneBack_f_s_shuff(i) = B_0_s_norm(TGM_test_shuff(i,2),TGM_test_shuff(i,3));%/possible_flight_num;
end


TwoBack_f = [];
for i=1:length(TGM_test)
    temp = ((TGM_test(i,2)-1)*li) + TGM_test(i,1);
    TwoBack_f(i) = B_1_norm(temp,TGM_test(i,3));
end

TwoBack_f_shuff = [];
for i=1:length(TGM_test_shuff)
    temp = ((TGM_test_shuff(i,2)-1)*li) + TGM_test_shuff(i,1);
    TwoBack_f_shuff(i) = B_1_norm(temp,TGM_test_shuff(i,3));
end

TwoBack_f_s = [];
for i=1:length(TGM_test)
    temp = ((TGM_test(i,2)-1)*li) + TGM_test(i,1);
    TwoBack_f_s(i) = B_1_s_norm(temp,TGM_test(i,3));
end


TwoBack_f_s_shuff = [];
for i=1:length(TGM_test_shuff)
    temp = ((TGM_test_shuff(i,2)-1)*li) + TGM_test_shuff(i,1);
    TwoBack_f_s_shuff(i) = B_1_s_norm(temp,TGM_test_shuff(i,3));
end

%%  KS test for OneBack versus Twoback subsetted and shuffled
KS_sample = 20;

clear h_s p_s HS_2 PS_2 HS_1 PS_1 TwoBackMatrix
TwoBack_Matrix = [TwoBack_f_s',TwoBack_f_s_shuff'];
for i=2:floor(length(TGM_test_shuff)/KS_sample)
    %indexes1 = randsample([1:size(TwoBack_f,2)],30);
    indexes1 = [(i-1)*KS_sample:(i-1)*KS_sample+KS_sample];
    [h_s_1,p_s_1] = kstest2(TwoBack_Matrix(indexes1,1),TwoBack_Matrix(indexes1,2),'Alpha',0.01/1000);
    HS_1(i) = h_s_1;
    PS_1(i) = p_s_1;
end
disp(strcat(num2str(sum(HS_1))," "," of 20 kstests of TwoBack_f vs. TwoBack_f_shuff were significantly different"));

clear h_s p_s HS_2 PS_2 HS_1 PS_1 OneBackMatrix
OneBack_Matrix = [OneBack_f_s',OneBack_f_s_shuff'];
for i=2:floor(length(TGM_test_shuff)/KS_sample)
    %indexes1 = randsample([1:size(OneBack_f,2)],30);
    indexes1 = [(i-1)*KS_sample:(i-1)*KS_sample+KS_sample];
    [h_s_1,p_s_1] = kstest2(OneBack_Matrix(indexes1,1),OneBack_Matrix(indexes1,2),'Alpha',0.01/1000);
    HS_1(i) = h_s_1;
    PS_1(i) = p_s_1;
end
disp(strcat(num2str(sum(HS_1))," "," of 20 kstests of OneBack_f vs. OneBack_f_shuff were significantly different"));

% OneBack_f vs. TwoBack_f
clear h_s p_s HS_2 PS_2 HS_1 PS_1 XBack_Matrix
XBack_Matrix = [OneBack_f_s',TwoBack_f_s'];
for i=2:floor(length(TGM_test_shuff)/KS_sample)
    %indexes1 = randsample([1:size(XBack_Matrix,1)],30);
    indexes1 = [(i-1)*KS_sample:(i-1)*KS_sample+KS_sample];
    [h_s_1,p_s_1] = kstest2(XBack_Matrix(indexes1,1),XBack_Matrix(indexes1,2),'Alpha',0.01/100);
    %figure(); hold on; plot(XBack_Matrix(indexes1,1)); plot(XBack_Matrix(indexes1,2));
    HS_1(i) = h_s_1;
    PS_1(i) = p_s_1;
end
disp(strcat(num2str(sum(HS_1))," "," of 20 kstests of OneBack_f vs. TwoBack_f were significantly different"));

% OneBack_f vs. TwoBack_f_shuff
clear h_s p_s HS_2 PS_2 HS_1 PS_1 XBack_Matrix
XBack_Matrix = [OneBack_f_s',TwoBack_f_s_shuff'];
for i=2:floor(length(TGM_test_shuff)/KS_sample)
    %indexes1 = randsample([1:size(XBack_Matrix,1)],30);
    indexes1 = [(i-1)*KS_sample:(i-1)*KS_sample+KS_sample];
    [h_s_1,p_s_1] = kstest2(XBack_Matrix(indexes1,1),XBack_Matrix(indexes1,2),'Alpha',0.01/100);
    HS_1(i) = h_s_1;
    PS_1(i) = p_s_1;
end
disp(strcat(num2str(sum(HS_1))," "," of 20 kstests of OneBack_f vs. TwoBack_f_shuff were significantly different"));
figure(); hold on; 
plot(cumsum(OneBack_f_s)); plot(cumsum(TwoBack_f_s_shuff));
