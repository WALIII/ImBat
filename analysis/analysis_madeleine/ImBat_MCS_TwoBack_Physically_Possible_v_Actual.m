function [h,p] = ImBat_MCS_TwoBack_Physically_Possible_v_Actual(c_s_34,pruned_dataset,takeoff_locations_cpu,outliers,varargin)

% Function to compare the conditional probability of the next flight conditioned on the takeoff location
% with the conditional probability of the next flight conditioned on the takeoff location and previous flight type. 
% For now there are two options for running this function:
%   1. Use the by-eye takeoff location vector (TPosition_vec, below)
%   2. Use the by-code takeoff location vector (Takeoff_Position_CPU)

% defaults
simulated_data_check = 1;
EC1=0;

% load in real arguments
while ~isempty(varargin)
    switch lower(varargin{1})
    	case 'simulated_data_check'
        	simulated_data_check = varargin{2};
        case 'exclude_cluster_1'
            EC1 = varargin{2};
        otherwise
        	error(['Unexpected option: ' varargin{1}])
    end
	varargin(1:2) = [];
end

% Make a 1 x num_clusters vector.
Cluster_vec = [1:size(unique(c_s_34),1)];
TPosition_vec = takeoff_locations_cpu;
  
DS = c_s_34(pruned_dataset);
DS_og = c_s_34;

% Determine which preceeding flights are incorrect due to pruning
pruned_index = {};
for i=2:size(c_s_34,1)
    pruned_set = size(outliers(outliers < i),2);
    pruned_index{i} = i-pruned_set;
end

has_incorrect_previous_flight = []; pairs_of_incorrect_previous_flights = [];
pruned_index_array = cell2mat(pruned_index);
for m=1:size(pruned_index_array,2)-2    
    if m==1 & ismember(1,outliers)
        continue;
    elseif pruned_index_array(m) == pruned_index_array(m+1)
        has_incorrect_previous_flight = [has_incorrect_previous_flight,pruned_index_array(m+1)];
        pairs_of_incorrect_previous_flights = [pairs_of_incorrect_previous_flights;[DS(pruned_index_array(m+1)-1),DS(pruned_index_array(m+1))]];
    end
end

% Find the takeoff position vector the same length as c_s_34 vector
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
            %landing_position(i) = LPosition_vec(j);
            takeoff_position(i) = TPosition_vec(j);
        end
    end
end

% For every unique takeoff position, get how many of that takeoff position
% there are (the prior distribution of takeoff positions)
for i=1:size(unique(takeoff_position(~isnan(takeoff_position))),2)
    %LPosition_P(i) = size(landing_position(landing_position==i),2);
    TPosition_P(i) = size(takeoff_position(takeoff_position==i),2);
end
% Get the probability of each takeoff position
Takeoff_Position_Prob = TPosition_P/size(takeoff_position,2);
if strcmp(EC1,'1')
    TPosition_P = TPosition_P(2:end);
    Takeoff_Position_Prob = TPosition_P/size(takeoff_position(takeoff_position~=1),2);  
end

%% Create A, which is unique positions by unique flight types, to store the
% number of times a flight type and a takeoff position occur together. This
% matrix will have many zero entries, do i keep these?... Ask Julie.

li = size(unique(takeoff_position(~isnan(takeoff_position))),2);
lk = size(unique(DS),1);
lxx = size(unique(takeoff_position(~isnan(takeoff_position))),2);

clear A;
A = zeros(li,lk);
for i=1:li
    for m=1:size(DS,1)
        if takeoff_position(m) == i
            for j=1:lk
                if DS(m) == j
                    A(i,j) = A(i,j)+1;
                end
            end
        end
    end
end

% Calculate the probability of a given flight, conditioning on takeoff
% location. All rows sum to 1.
AA = [];
for i=1:size(A,1)
    row_of_concern = A(i,:);
    for j=1:size(row_of_concern,2)
        AA(i,j) = A(i,j)./sum(row_of_concern);
    end
end

% If we exclude cluster 1 (exclude column 1), renormalize.
if strcmp(EC1,'1')
    A = A(:,2:end);
    AA = [];
    for i=1:size(A,1)
        row_of_concern = A(i,:);
        for j=1:size(row_of_concern,2)
            AA(i,j) = A(i,j)./sum(row_of_concern);
        end
    end
end

AA(isnan(AA))=0;

%% Now construct B, conditioning also on the previous flight.
clear B;
B = zeros(li*lk,lk);
tt=0;
for i=1:li
    for k=1:lk
        for m=2:size(DS,1)   
            if takeoff_position(m) == i & DS(m-1) == k & ~ismember(m,has_incorrect_previous_flight)
                for j=1:size(unique(DS),1)
                    if DS(m) == j
                        %B(25*(i-1)+k,j) =B(25*(i-1)+k,j)+1;
                        B(size(unique(DS),1)*(i-1)+k,j) =B(size(unique(DS),1)*(i-1)+k,j)+1;
                    end
                end
            elseif takeoff_position(m) == i & DS(m-1) == k & ismember(m,has_incorrect_previous_flight)
                tt = tt+1;
                for j=1:size(unique(DS),1)
                    if c_s_34(outliers(tt)) == j
                        B(size(unique(DS),1)*(i-1)+k,j) = B(size(unique(DS),1)*(i-1)+k,j)+1;
                    end
                end
            end
        end
    end
end


% Calculate the conditional probability; each row sums to 1
BB = zeros(size(B,1),size(B,2));
for i=1:size(B,1)
    row_of_concern = B(i,:);
    for j=1:size(row_of_concern,2)
        BB(i,j) = B(i,j)/sum(row_of_concern);
    end
end

BB(isnan(BB))=0;

%% Now construct CC, conditioning also on the previous flight and previous takeoff location
clear C;

C = zeros(li*lk*lxx,lk);
tt=0;
for i=1:li % Takeoff pos
    for k=1:lk                                      % Prev flight type 
        for xx = 1:lxx+1 % Prev Takeoff pos
            disp(((i-1)*li) + ((k-1)*lk) + xx)
            for m=2:size(DS,1)
                if takeoff_position(m) == i & DS(m-1) == k & takeoff_position(m-1) == xx & ~ismember(m,has_incorrect_previous_flight)
                    for j=1:lk
                        if DS(m) == j
                            C(((i-1)*li*lk) + ((k-1)*lxx) + xx, j) = C(((i-1)*li*lk) + ((k-1)*lxx) + xx, j)+ 1;
                        end
                    end
                elseif takeoff_position(m) == i & DS(m-1) == k & takeoff_position(m-1) == xx & ismember(m,has_incorrect_previous_flight)
                    tt = tt+1;
                    for j=1:size(unique(DS),1)
                        if c_s_34(outliers(tt)) == j
                            C(((i-1)*li*lk) + ((k-1)*lxx) + xx, j) = C(((i-1)*li*lk) + ((k-1)*lxx) + xx, j)+ 1;
                        end
                    end
                end
            end
        end
    end
end


% Calculate the conditional probability; each row sums to 1
CC = zeros(size(C,1),size(C,2));
for i=1:size(C,1)
    row_of_concern = C(i,:);
    for j=1:size(row_of_concern,2)
        CC(i,j) = C(i,j)/sum(row_of_concern);
    end
end

CC(isnan(CC))=0;

%% 
CTR = [];
for i=1:size(unique(takeoff_position(~isnan(takeoff_position))),2)+1
    to_add = (i-1)*size(unique(DS),1);
    CTR = [CTR, to_add];
end
CTR = CTR+1;
CTR = CTR(1:end-1);

% If we exclude cluster 1 (exclude column 1), renormalize.
if strcmp(EC1,'1')
    B(CTR,:) = [];
    BB = [];
    for i=1:size(B,1)
        row_of_concern = B(i,:);
        for j=1:size(row_of_concern,2)
            BB(i,j) = B(i,j)./sum(row_of_concern);
        end
    end
end

BB(isnan(BB))=0;

%% Calculate the prior distribution of positions

clear TPPC;
% Calculate the prior distribution of positions/previous flight types   
TPPC = zeros(li*lk,1);
for i=1:li
    for k=1:lk
        for m=2:size(DS,1)
            if takeoff_position(m) == i & DS(m-1) == k
                %TPPC(max(DS)*(i-1)+k) = TPPC(max(DS)*(i-1)+k)+1;
                TPPC(size(unique((DS)),1)*(i-1)+k) = TPPC(size(unique((DS)),1)*(i-1)+k)+1;
            end
        end
    end
end
                
Takeoff_Position_Cond_Prob_B = TPPC./sum(TPPC);
% Remove rows with cluster 1 if needed
if strcmp(EC1,'1')
    TPPC(CTR,:) = [];
    Takeoff_Position_Cond_Prob_B = TPPC./sum(TPPC);
end

clear TPPC;
% Calculate the prior distribution of previous flight type and positions 
TPPC = zeros(li*lk*lxx,1);
for i=1:li
    for k=1:lk
        for xx = 1:lxx
            for m=2:size(DS,1)
                if takeoff_position(m) == i & DS(m-1) == k & takeoff_position(m-1) == xx
                    TPPC(((i-1)*li*lk) + ((k-1)*lxx) + xx) = TPPC(((i-1)*li*lk) + ((k-1)*lxx) + xx)+1;
                end
            end
        end
    end
end
               
Takeoff_Position_Cond_Prob_C = TPPC./sum(TPPC);

% Remove rows with cluster 1 if needed
if strcmp(EC1,'1')
    TPPC(CTR,:) = [];
    Takeoff_Position_Cond_Prob_C = TPPC./sum(TPPC);
end

%% Joints
% AA
Joint_AA = zeros(size(A,1),size(A,2));

for i=1:size(Takeoff_Position_Prob,2)
    t_p_i = Takeoff_Position_Prob(i);
    AA_row = AA(i,:);
    Joint_AA(i,:) = AA_row'*t_p_i;
end

Joint_AA(isnan(Joint_AA))=0;
Joint_AA_vec = reshape(Joint_AA,1,size(Joint_AA,1)*size(Joint_AA,2));

% BB
%Joint_BB = zeros(size(unique((DS)),1)*size(unique((DS)),1),size(unique((DS)),1));
%Joint_BB = zeros(size(unique(takeoff_position(~isnan(takeoff_position))),2)*size(unique((DS)),1),size(unique((DS)),1));
Joint_BB = zeros(size(BB,1),size(BB,2));
for i=1:size(Takeoff_Position_Cond_Prob_B,1)
    t_p_i = Takeoff_Position_Cond_Prob_B(i);
    BB_row = BB(i,:);
    Joint_BB(i,:) = BB_row'*t_p_i;
end

Joint_BB_vec = reshape(Joint_BB,1,size(Joint_BB,1)*size(Joint_BB,2));

% CC
Joint_CC = zeros(size(CC,1),size(CC,2));
for i=1:size(Takeoff_Position_Cond_Prob_C,1)
    t_p_i = Takeoff_Position_Cond_Prob_C(i);
    CC_row = CC(i,:);
    Joint_CC(i,:) = CC_row'*t_p_i;
end

Joint_CC_vec = reshape(Joint_CC,1,size(Joint_CC,1)*size(Joint_CC,2));

%% KS Test
clear p h;
[h,p] = kstest2(Joint_BB_vec,Joint_CC_vec); 

end