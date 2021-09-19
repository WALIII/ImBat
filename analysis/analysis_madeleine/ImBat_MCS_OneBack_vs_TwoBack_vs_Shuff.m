% Try
% 1. shuffle
% 2. test
% 3. repeat

[TGM_sorted,TGM_ss] = sort(TGM(:,2));
temp = TGM(:,1);
TGM_sorted_1 = temp(TGM_ss);
temp = TGM(:,3);
TGM_sorted_2 = temp(TGM_ss);
TGM_sorted_full = [TGM_sorted_1,TGM_sorted,TGM_sorted_2];
HS_1 = [];
HS_1_shuff = [];
co=19; %li;
for ii = 1:1000
    TGM_sorted_full_shuff = [];
    for i=1:length(subsegment_matrix)    
        temp = subsegment_matrix{i};
        temp_1 = temp(:,1);
        temp_1_shuff = temp_1(randperm(length(temp_1)));
        temp_shuff = [temp_1_shuff,temp(:,2),temp(:,3)];
        TGM_sorted_full_shuff = [TGM_sorted_full_shuff;temp_shuff];
    end
    
    f=3;  
%     TGM_test = []; TGM_test_shuff = [];
%     for i=1:length(TGM_sorted_full)
%         if TGM_sorted_full(i,2) == f
%             TGM_test = [TGM_test;TGM_sorted_full(i,:)];
%         end
%         if TGM_sorted_full_shuff(i,2) == f
%             TGM_test_shuff = [TGM_test_shuff;TGM_sorted_full_shuff(i,:)];
%         end
%     end
    TGM_test = TGM_sorted_full((TGM_sorted_full(:,2)==f),:);
    TGM_test_shuff = TGM_sorted_full_shuff((TGM_sorted_full_shuff(:,2)==f),:);

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

    TwoBack_f = [];
    for i=1:length(TGM_test)
        temp = ((TGM_test(i,2)-1)*(co+1)) + TGM_test(i,1);
        TwoBack_f(i) = B_1_norm(temp,TGM_test(i,3));
    end

    TwoBack_f_shuff = [];
    for i=1:length(TGM_test_shuff)
        temp = ((TGM_test_shuff(i,2)-1)*(co+1)) + TGM_test_shuff(i,1);
        TwoBack_f_shuff(i) = B_1_norm(temp,TGM_test_shuff(i,3));
    end


    %%  KS test for OneBack versus Twoback subsetted and shuffled
    KS_sample = 20;

    % OneBack_f vs. TwoBack_f
    clear h_s_1 XBack_Matrix indexes1
    XBack_Matrix = [OneBack_f',TwoBack_f'];
    %for i=2:floor(length(TGM_test_shuff)/KS_sample)
    indexes1 = randsample([1:size(XBack_Matrix,1)],15);
    %figure(); hold on; title('OneBack vs. TwoBack'); plot(XBack_Matrix(indexes1,1)); plot(XBack_Matrix(indexes1,2));
    [h_s_1,p_s_1] = kstest2(XBack_Matrix(indexes1,1),XBack_Matrix(indexes1,2),'Alpha',0.01/1000);
    %[h_s_1,p_s_1] = kstest2(XBack_Matrix(:,1),XBack_Matrix(:,2),'Alpha',0.01/1000);
    HS_1(ii) = h_s_1;

    % OneBack_f vs. TwoBack_f_shuff
    clear h_s_1 XBack_Matrix indexes1
    XBack_Matrix = [OneBack_f',TwoBack_f_shuff'];
    indexes1 = randsample([1:size(XBack_Matrix,1)],15);
    %figure(); hold on; title('OneBack vs. Shuffle'); plot(XBack_Matrix(indexes1,1)); plot(XBack_Matrix(indexes1,2));
    [h_s_1,p_s_1] = kstest2(XBack_Matrix(indexes1,1),XBack_Matrix(indexes1,2),'Alpha',0.01/1000);
    %[h_s_1,p_s_1] = kstest2(XBack_Matrix(:,1),XBack_Matrix(:,2),'Alpha',0.01/1000);
    HS_1_shuff(ii) = h_s_1;
    
    %figure(); hold on; plot(cumsum(OneBack_f),'r'); plot(cumsum(TwoBack_f),'g'); plot(cumsum(TwoBack_f_shuff),'b');

    
end

figure(); hold on; plot(cumsum(OneBack_f),'r'); plot(cumsum(TwoBack_f),'g'); plot(cumsum(TwoBack_f_shuff),'b');
disp(sum(HS_1))
disp(sum(HS_1_shuff))
