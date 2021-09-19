% Script to test if: 
% P(F(t+1), F(t)) == P(F(t+1), F(t), F(t-1))

clear all; close all;
load('c_s_34.mat');

KS_sample = 20;
% Flight cluster:
f = 5;

% Use all data if not doing the rescaling above
%c_s_34_seg = subset_c_s_34%c_s_34;%(1000:end-1000);
%c_s_34_seg_r = subset_c_s_34%c_s_34;%(1000:end-1000);

% Use subset of data from 57 to 67
c_s_34_seg = c_s_34;
c_s_34_seg_r = c_s_34;


    
% Plot prior distribution of flightpaths
figure(); hold on; histogram(c_s_34_seg_r); xlabel("Flight Type"); ylabel("Frequency"); title("Prior Distribution of Flight Types");

% Extract every 3-gram possible from the data using sliding window
TGM_all = [];
window_size=2;
for i=1:length(c_s_34_seg_r)-window_size
    three_gram = c_s_34_seg_r(i:i+window_size);
    TGM_all(i,:) = three_gram;
end

% Calculate co based on which flight types make up >5% of the data
five_percent = round(size(c_s_34,1)*(0.05));
length(c_s_34(c_s_34 > 50));
% Remove all 3-grams with less than co in them
co=16;%25;
li = size(unique(TGM_all),1);
TGM=[];
for i=1:size(TGM_all,1)
    row = TGM_all(i,:);
    if (sum(row > co+1) >= 1) | ismember(1,row)
        continue;
    else
        TGM = [TGM;row];
    end
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

% Calculate normalization factor (how many different flight types occur
% before a flight type type f.
possible_flight_num_temp = TGM(TGM(:,2)==f,:);
possible_flight_num = length(unique(possible_flight_num_temp(:,1)));

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
for i=1:length(TGM)
    OneBack_s(i) = B_0_s_norm(TGM(i,2),TGM(i,3));
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
for i=1:length(TGM)
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
B_0_trimmed = B_0(2:co+1,2:co+1);
B_0_norm = [];
for i=1:size(B_0_trimmed,1)
    row = B_0_trimmed(i,:);
    row_scaled = row./possible_flight_num;
    row_sum = row./sum(row);
    B_0_norm(i,:) = row_sum;
end
B_0_norm(isnan(B_0_norm))=0;

OneBack = [];
for i=1:length(TGM)
    OneBack(i) = B_0_norm(TGM(i,2)-1,TGM(i,3)-1);
end

B_1 = zeros(size(B_0_norm,1)*size(B_0_norm,1),size(B_0_norm,1));
for i=2:(co+1)   
    for j=2:(co+1) 
        for m=3:length(c_s_34_seg_r)                          
            for k=2:(co+1)                                  
                if c_s_34_seg_r(m-2) == j & c_s_34_seg_r(m-1)== i & c_s_34_seg_r(m) == k    
                    B_1((i-2)*(co) + (j-1),k-1) = B_1((i-2)*(co) + (j-1),k-1) + 1;                 
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
for i=1:length(TGM)
    P_B_1_idx = ((TGM(i,2)-2)*co) + (TGM(i,1)-1);
    TwoBack(i) = B_1_norm(P_B_1_idx,TGM(i,3)-1);
end

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
for i=2:length(temp_unique)+1
    subsegs = [];
    for j=1:length(TGM_sorted_full)
        if TGM_sorted_full(j,2) == temp_unique(i-1)
            subsegs = [subsegs; TGM_sorted_full(j,:)];
        end
    end
    subsegment_matrix{i-1} = subsegs;
end

TGM_sorted_full_shuff = [];
for i=1:length(subsegment_matrix)    
    temp = subsegment_matrix{i};
    temp_1 = temp(:,1);
    temp_1_shuff = temp_1(randperm(length(temp_1)));
    temp_shuff = [temp_1_shuff,temp(:,2),temp(:,3)];
    TGM_sorted_full_shuff = [TGM_sorted_full_shuff;temp_shuff];
end

%% Now do the test for the shuffled data matrix
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
    OneBack_f(i) = B_0_norm(TGM_test(i,2)-1,TGM_test(i,3)-1);
end

OneBack_f_shuff = [];
for i=1:length(TGM_test_shuff)
    OneBack_f_shuff(i) = B_0_norm(TGM_test_shuff(i,2)-1,TGM_test_shuff(i,3)-1);
end

% OneBack_f_s = [];
% for i=1:length(TGM_test_shuff)
%     OneBack_f_s(i) = B_0_s_norm(TGM_test(i,2)-1,TGM_test(i,3)-1);
% end
% 
% OneBack_f_s_shuff = [];
% for i=1:length(TGM_test_shuff)
%     OneBack_f_s_shuff(i) = B_0_s_norm(TGM_test_shuff(i,2)-1,TGM_test_shuff(i,3)-1);
% end

TwoBack_f = [];
for i=1:length(TGM_test)
    temp = ((TGM_test(i,2)-2)*co) + (TGM_test(i,1)-1);
    TwoBack_f(i) = B_1_norm(temp,TGM_test(i,3)-1);
end

TwoBack_f_shuff = [];
for i=1:length(TGM_test_shuff)
    temp = ((TGM_test_shuff(i,2)-2)*co) + (TGM_test_shuff(i,1)-1);
    TwoBack_f_shuff(i) = B_1_norm(temp,TGM_test_shuff(i,3)-1);
end

% TwoBack_f_s = [];
% for i=1:length(TGM_test)
%     temp = ((TGM_test(i,2)-1)*20) + TGM_test(i,1);
%     TwoBack_f_s(i) = B_1_s_norm(temp,TGM_test(i,3)-1);
% end
% 
% TwoBack_f_s_shuff = [];
% for i=1:length(TGM_test_shuff)
%     temp = ((TGM_test_shuff(i,2)-1)*li) + TGM_test_shuff(i,1);
%     TwoBack_f_s_shuff(i) = B_1_s_norm(temp,TGM_test_shuff(i,3)-1);
% end

%%  KS test for OneBack versus Twoback subsetted and shuffled

% OneBack_f vs. TwoBack_f
clear h_s p_s HS_2 PS_2 HS_1 PS_1 XBack_Matrix
XBack_Matrix = [OneBack_f',TwoBack_f'];
test_num = floor(length(TGM_test_shuff)/KS_sample);
test_num=5000;
for i=2:test_num
    indexes1 = randsample([1:size(XBack_Matrix,1)],KS_sample);
    %indexes1 = [(i-1)*KS_sample:(i-1)*KS_sample+KS_sample];
    [h_s_1,p_s_1] = kstest2(XBack_Matrix(indexes1,1),XBack_Matrix(indexes1,2),'Alpha',0.05/test_num);
    %figure(); hold on; plot(XBack_Matrix(indexes1,1)); plot(XBack_Matrix(indexes1,2));
    HS_1(i) = h_s_1;
    PS_1(i) = p_s_1;
end
disp(strcat(num2str(sum(HS_1))," "," of "," ",num2str(test_num)," ","kstests of OneBack_f vs. TwoBack_f were significantly different"));
figure(); hold on; 
plot(cumsum(OneBack_f)); plot(cumsum(TwoBack_f)); title(strcat("OneBack_f v.s. TwoBack_f where"," ", num2str(f)," ","is flight"));
figure(); hold on; histogram(PS_1); title("OneBack v.s. Twoback")

% OneBack_f vs. TwoBack_f_shuff
clear h_s p_s HS_2 PS_2 HS_1 PS_1 XBack_Matrix
XBack_Matrix = [OneBack_f',TwoBack_f_shuff'];
for i=2:test_num
    indexes1 = randsample([1:size(XBack_Matrix,1)],KS_sample);
    %indexes1 = [(i-1)*KS_sample:(i-1)*KS_sample+KS_sample];
    [h_s_1,p_s_1] = kstest2(XBack_Matrix(indexes1,1),XBack_Matrix(indexes1,2),'Alpha',0.05/test_num);
    HS_1(i) = h_s_1;
    PS_1(i) = p_s_1;
end
disp(strcat(num2str(sum(HS_1))," "," of "," ",num2str(test_num)," ","kstests of OneBack_f vs. TwoBack_f shuffled were significantly different"));
figure(); hold on; 
plot(cumsum(OneBack_f)); plot(cumsum(TwoBack_f_shuff)); title(strcat("OneBack_f v.s. TwoBack_f shuffled where"," ", num2str(f)," ","is flight"));
figure(); hold on; histogram(PS_1); title("OneBack v.s. Twoback Shuff")

% for xxx = 1:25
%     KS_sample = xxx+10;
%     % OneBack_f vs. TwoBack_f
%     clear h_s p_s HS_2 PS_2 HS_1 PS_1 XBack_Matrix
%     XBack_Matrix = [OneBack_f',TwoBack_f'];
%     test_num = floor(length(TGM_test_shuff)/KS_sample);
%     test_num=1000;
%     for i=2:test_num
%         indexes1 = randsample([1:size(XBack_Matrix,1)],KS_sample);
%         %indexes1 = [(i-1)*KS_sample:(i-1)*KS_sample+KS_sample];
%         [h_s_1,p_s_1] = kstest2(XBack_Matrix(indexes1,1),XBack_Matrix(indexes1,2),'Alpha',0.05/test_num);
%         %figure(); hold on; plot(XBack_Matrix(indexes1,1)); plot(XBack_Matrix(indexes1,2));
%         HS_1(i) = h_s_1;
%         PS_1(i) = p_s_1;
%     end
%     cumulative_HS_1_1(xxx) = sum(HS_1)/test_num;
%     mean_PS_1_1(xxx) = mean(PS_1);
%     
%     % OneBack_f vs. TwoBack_f_shuff
%     clear h_s p_s HS_2 PS_2 HS_1 PS_1 XBack_Matrix
%     XBack_Matrix = [OneBack_f',TwoBack_f_shuff'];
%     for i=2:test_num
%         indexes1 = randsample([1:size(XBack_Matrix,1)],KS_sample);
%         %indexes1 = [(i-1)*KS_sample:(i-1)*KS_sample+KS_sample];
%         [h_s_1,p_s_1] = kstest2(XBack_Matrix(indexes1,1),XBack_Matrix(indexes1,2),'Alpha',0.05/test_num);
%         HS_1(i) = h_s_1;
%         PS_1(i) = p_s_1;
%     end
%     cumulative_HS_1_2(xxx) = sum(HS_1)/test_num;
%     mean_PS_1_2(xxx) = mean(PS_1);
% end
% figure(); hold on; plot(cumulative_HS_1_1); plot(cumulative_HS_1_2); yline(0.05);
% xlabel("KS Test sample size");
% ylabel("Significance of Null and Alt");
% title("H values as a function of KS Test sample sizes");
% 
% figure(); hold on; plot(mean_PS_1_1); plot(mean_PS_1_2); yline(0.05);
% xlabel("KS Test sample size");
% ylabel("Significance of Null and Alt");
% title("P values as a function of KS Test sample sizes");

    