function [FnormM,FnormD] = ImBat_MCS_Calculate_Simulated_T_Matrix(S,num_clust,day_select)

    %% T matrix for simulated sequences

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

        clear F;
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
                 %   fp = sum(bpc==j);
                 %   FP(i,j) = fp;
                %else
                %    FP(i,j) = 0;
                %end
            end
        end

        % Fnorm is the transition probability of x to y
        % FPnorm is the transition probability from y to x

        % Construct naive probability matrix
        Fnorm34sim = []; %FPnorm34sim = [];
        for i=1:num_clust
            row = F(i,:);
            %rowP = FP(i,:);
            Fnorm34sim(i,:) = row./sum(row);
            %FPnorm34sim(i,:) = rowP./sum(rowP);
        end

        %figure(); h = heatmap(Fnorm34, 'Colormap', copper); h.title(strcat("Simulated Sequence #",num2str(m),"; Day ", num2str(day_select)));
        FnormD(:,:,m) = Fnorm34sim;
    end

    FnormM = nanmean(FnormD,3);
    figure(); h = heatmap(FnormM, 'Colormap', copper); h.title(strcat("Mean Simulated Sequence; Day ", num2str(day_select)));

end    