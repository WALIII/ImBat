function [H_SIM,P_SIM] = ImBat_MCS_Physically_Possible_v_Actual_simulated(S,pruned_dataset,takeoff_locations_cpu,varargin)

    % Defaults
    EC1 = 0;

    % Do this for simulated data to check
    for x=1:99
        clear takeoff_position LPosition_P TPosition_P TPosition_vec TPPC;
        Sdata = S(x,:)';
        DS = Sdata;
        
        % Vector 1;num_clusters
        Cluster_vec = [unique(Sdata)];
        TPosition_vec = takeoff_locations_cpu;
        
        clear landing_position takeoff_position;
        takeoff_position = zeros(1,size(Sdata,1));
        if length(Cluster_vec) > length(TPosition_vec)
            smaller_size = length(TPosition_vec);
        else
            smaller_size = length(Cluster_vec);
        end
        for j=1:smaller_size
            for i=1:size(Sdata,1)
                if Sdata(i) == Cluster_vec(j)
                    %landing_position(i) = LPosition_vec(j);
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
        %Landing_Position_Prob = LPosition_P/size(landing_position,2);
        Takeoff_Position_Prob = TPosition_P/size(takeoff_position,2);
        if strcmp(EC1,'1')
            TPosition_P = TPosition_P(2:end)
            Takeoff_Position_Prob = TPosition_P/size(takeoff_position(takeoff_position~=1),2);  
        end


        % Create A, which is unique positions by unique flight types, to store the
        % number of times a flight type and a takeoff position occur together. This
        % matrix will have many zero entries, do i keep these?... Ask Julie.
        clear A;
        A = zeros(size(unique(takeoff_position(~isnan(takeoff_position))),2),size(unique(DS),1));
        for i=1:size(unique(takeoff_position(~isnan(takeoff_position))),2)
            for m=1:size(DS,1)
                if takeoff_position(m) == i
                    for j=1:size(unique(DS),1)
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

        % Now construct AB but conditioning also on the previous flight.
        B = zeros(size(unique(takeoff_position(~isnan(takeoff_position))),2)*size(unique(DS),1),size(unique(DS),1));
        for i=1:size(unique(takeoff_position(~isnan(takeoff_position))),2)
            for k=1:size(unique(DS),1)
                for m=2:size(DS,1)
                    if takeoff_position(m) == i & DS(m-1) == k
                        for j=1:size(unique(DS),1)
                            if DS(m) == j
                                %B(25*(i-1)+k,j) =B(25*(i-1)+k,j)+1;
                                B(size(unique(DS),1)*(i-1)+k,j) =B(size(unique(DS),1)*(i-1)+k,j)+1;
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
        TPPC = zeros(size(unique(takeoff_position(~isnan(takeoff_position))),2)*size(unique(DS),1),1);
        for i=1:size(unique(takeoff_position(~isnan(takeoff_position))),2)
            for k=1:size(unique(DS),1)
                for m=2:size(DS,1)
                    if takeoff_position(m) == i & DS(m-1) == k
                        %TPPC(max(DS)*(i-1)+k) = TPPC(max(DS)*(i-1)+k)+1;
                        TPPC(size(unique((DS)),1)*(i-1)+k) = TPPC(size(unique((DS)),1)*(i-1)+k)+1;
                    end
                end
            end
        end

        Takeoff_Position_Cond_Prob = TPPC./sum(TPPC);
        % Remove rows with cluster 1 if needed
        if strcmp(EC1,'1')
            TPPC(CTR,:) = [];
            Takeoff_Position_Cond_Prob = TPPC./sum(TPPC);
        end

        %Joints
        % AA
        Joint_AA = zeros(size(A,1),size(A,2));
        %Joint_AA = zeros(size(unique(takeoff_position(~isnan(takeoff_position))),2),size(unique((DS)),1));
        for i=1:size(Takeoff_Position_Prob,2)
            t_p_i = Takeoff_Position_Prob(i);
            AA_row = AA(i,:);
            Joint_AA(i,:) = AA_row'*t_p_i;
        end

        Joint_AA(isnan(Joint_AA))=0;
        Joint_AA_vec = reshape(Joint_AA,1,size(Joint_AA,1)*size(Joint_AA,2));

        % BB
        Joint_BB = zeros(size(B,1),size(B,2));
        for i=1:size(Takeoff_Position_Cond_Prob,1)
            t_p_i = Takeoff_Position_Cond_Prob(i);
            BB_row = BB(i,:);
            Joint_BB(i,:) = BB_row'*t_p_i;
        end

        Joint_BB_vec = reshape(Joint_BB,1,size(Joint_BB,1)*size(Joint_BB,2));

        % KS Test
        [h_sim,p_sim] = kstest2(Joint_AA_vec,Joint_BB_vec);
        H_sim(x) = h_sim;
        P_sim(x) = p_sim;
    end
    H_SIM = H_sim;
    P_SIM = P_sim;
end