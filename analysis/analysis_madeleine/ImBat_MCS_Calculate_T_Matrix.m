function [Fnorm34,c_s_Tnorm_34,c_s_T_34,OG_Fnorm] = ImBat_MCS_Calculate_T_Matrix(c_s_34,num_clust)

    % Get the T matrix for the Markov Model
    % Each flight type is a state

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

    % Save this Fnorm 
    OG_Fnorm = Fnorm34;

    Fnorm34_first5 = Fnorm34(1:5,1:5);
    Fnorm34_noone = Fnorm34(2:10,2:10);
    Fnorm34_first10 = Fnorm34(1:10,1:10);
    FPnorm34_first10 = FPnorm34(1:10,1:10);

    figure(); h = heatmap(Fnorm34, 'Colormap', copper); h.title("All Flight Cluster Transition Probability Matrix");
    figure(); h = heatmap(Fnorm34_first10, 'Colormap', copper);  h.title("First 10 Flight Cluster Transition Probability Matrix");

    %% Get the T matrix for the Markov Model
    % sStereotyped v Non-stereotyped is a state
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

    c_s_Tnorm_34(1,:) = c_s_T_34(1,:)./sum( c_s_T_34(1,:));
    c_s_Tnorm_34(2,:) = c_s_T_34(2,:)./sum( c_s_T_34(2,:));
    
    mc = dtmc(c_s_Tnorm_34); 
    mc.StateNames = [" " " "];
    figure(); 
    h = graphplot(mc);
    c = h.EdgeColor;
    h.EdgeColor = 'k';
    mcp = nonzeros(reshape(mc.P,[size(mc.P,1)*size(mc.P,1),1])).*8;
    h.LineWidth = mcp;
    h.MarkerSize = [(c_s_T_34(1,1)+c_s_T_34(1,2))/100,
                    (c_s_T_34(2,1)+c_s_T_34(2,2))/100];

end

