% Script to determine if 
% P(F(t+1), T(t)) == P(F(t+1), F(t), T(t))
% 
% color_scheme = [parula(10); parula(length(unique(c_s_34))-10)];
f = 3;
% L = [];
% for i=2:length(c_s_34)-1
%     if c_s_34(i) == f
%         L = [L;c_s_34(i-1),c_s_34(i),c_s_34(i+1),i];
%     end
% end
% 
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

%% Make takeoff_positions vector 
Cluster_vec = [1:size(unique(c_s_34),1)];
TPosition_vec = takeoff_locations_cpu;
DS = c_s_34;

clear landing_position takeoff_position LPosition_P TPosition_P;
takeoff_position = zeros(1,size(DS,1));
if length(Cluster_vec) > length(TPosition_vec)
    smaller_size = length(TPosition_vec);
else
    smaller_size = length(Cluster_vec);
end
for i=1:size(DS,1)
    for j=1:smaller_size
        if DS(i) == Cluster_vec(j)
            takeoff_position(i) = TPosition_vec(j);
        end
    end
end

% For every unique takeoff position, get how many of that takeoff position
% there are (the prior distribution of takeoff positions)
for i=1:size(unique(takeoff_position(~isnan(takeoff_position))),2)
    TPosition_P(i) = size(takeoff_position(takeoff_position==i),2);
end
% Get the probability of each takeoff position
Takeoff_Position_Prob = TPosition_P/size(takeoff_position,2);

li = length(unique(c_s_34));
lk = length(unique(takeoff_locations_cpu));

%% Calculate:
% P(F(t+1), T(t)) = P(F(t+1) | T(t)) * P(T(t))

% A = counts of P(F(t+1) | T(t))
% Row index is T(t), Column index is F(t+1)
A = zeros(lk,li);
for i=1:lk                                                  % Pick a takeoff location
    for m=1:length(c_s_34)                                  % For every Flight,
        for j=1:li                                          % For every flight type
            if takeoff_position(m)== i & c_s_34(m) == j    % If the takeoff location was of type i, and next flight was of type j, 
                A(i,j) = A(i,j)+1;                          % Add to the matrix
            end
        end
    end
end

% The total number of flights stored in A should equal the length of c_s_34-1
sum(sum(A)) == length(c_s_34)

% P_A = P(F(t+1) | T(t))
P_A = [];
for i=1:size(A,1)
    row = A(i,:);
    P_A(i,:) = row./sum(row);
end

P_A(isnan(P_A))=0;
sum(P_A,2);

% Counts of P(T(t))
A_prior = zeros(lk,1);
for i=1:lk                                          
    for m=1:length(c_s_34)                         
        if takeoff_position(m) == i    
            A_prior(i) = A_prior(i)+1;                  
        end
    end
end

% P_A = P(F(t))
P_A_prior = A_prior./(length(c_s_34));


% P(F(t+1), F(t)) = P(F(t+1) | F(t)) * P(F(t))
Joint_A = mtimes(P_A',P_A_prior);

%%
% B_1 = P(F(t+1) | T(t), F(t)) COUNTS
% Row index is (T(t), F(t)), Column index is F(t+1)
B_1 = zeros(lk*li,li);
for i=1:lk   
    for j=1:li 
        for m=2:length(c_s_34)                          
            for k=1:li                                  
                if c_s_34(m-1) == j & takeoff_position(m) == i & c_s_34(m) == k    
                    B_1((i-1)*li + j,k) = B_1((i-1)*li + j,k) + 1;                 
                end
            end
        end
    end
end

% The total number of flights stored in B should equal the length of c_s_34-2
sum(sum(B_1)) == length(c_s_34)-1;

% P_B_1 = P(F(t+1) | T(t), F(t))
P_B_1 = zeros(size(B_1,2),size(B_1,1));
for i=1:size(B_1,2)
    row = B_1(:,i);
    P_B_1(i,:) = row./sum(row);
end
P_B_1(isnan(P_B_1))=0;

sum(P_B_1,2);

% B_2 = P(F(t) | T(t)) COUNTS
B_2 = zeros(lk,li);
for i=1:lk   
    for m=2:length(c_s_34)                          
        for k=1:li                                  
            if c_s_34(m-1) == k & takeoff_position(m)== i  
                B_2(i,k) = B_2(i,k) + 1;                 
            end
        end
    end
end

% The total number of flights stored in B should equal the length of c_s_34-2
sum(sum(B_2)) == length(c_s_34)-1;

% P_B_2 = P(F(t) | T(t))
P_B_2 = zeros(size(B_2,1),size(B_2,2));
for i=1:size(B_2,1)
    row = B_2(i,:);
    P_B_2(i,:) = row./sum(row);
end
P_B_2(isnan(P_B_2))=0;

sum(P_B_2,2);

% B_3 = P(T(t)) COUNTS
B_3 = zeros(lk,1);
for i=1:lk                                          
    for m=1:length(c_s_34)                         
        if takeoff_position(m) == i    
            B_3(i) = B_3(i)+1;                  
        end
    end
end

% P_A = P(T(t))
P_B_3 = B_3./(length(c_s_34));

% P(F(t+1), T(t), F(t)) = P(F(t+1) | T(t)) , (F(t)) * P(T(t)) | (F(t)) * P(F(t))
Joint_temp = mtimes(P_B_2',P_B_3);
Joint_B = mtimes(P_B_1',Joint_temp);

%% KS Test 
[h,p] = kstest2(Joint_A,Joint_B);

%% P(Ft+1) | F(t) = f)
BB = [];
AA = P_A(f,:);
for i=1:size(P_A,1)
    BB(i,:) = P_B_1(:,(i-1)*size(P_A,1) + f)';
end

figure(); hold on; plot(BB); plot(AA,'LineWidth',2); title("Difference in probability of next flight, given that previous flight is 2, after conditioning on prev prev flight");














