function [legal_sequences_day34,weighted_legal_sequences_day34,WLS_day34,LS_day34,S] = ImBat_MCS_simulated_sequences(Fnorm34,c_s_34,weighted)

    % Calculate the T matrix for num_shuffle fake sequences.

    legal_sequences_day34 = cell(1,size(Fnorm34,1));
    weighted_legal_sequences_day34 = cell(1,size(Fnorm34,1));
    for i=1:size(Fnorm34,1)
        nonzero_idxs = find(Fnorm34(i,:) ~= 0);
        legal_sequences_day34{i} = nonzero_idxs;
        for j=1:size(nonzero_idxs,2)
            legal_sequences_day34{i} = [legal_sequences_day34{i},nonzero_idxs(j)];
        end
        legal_sequences_day34{i} = unique(legal_sequences_day34{i});
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
                if weighted == 1
                    tvec = WLS_day34{1};
                    seq_bloc = tvec(randi(length(tvec)));
                    new_seq = [new_seq;seq_bloc];
                else 
                    tvec = LS_day34{1};
                    seq_bloc = tvec(randi(length(tvec)));
                    new_seq = [new_seq;seq_bloc];
                end
            else
                if weighted == 1
                    tvec = WLS_day34{new_seq(j-1)};
                    seq_bloc = tvec(randi(length(tvec)));
                    new_seq = [new_seq;seq_bloc];
                else 
                    tvec = LS_day34{new_seq(j-1)};
                    seq_bloc = tvec(randi(length(tvec)));
                    new_seq = [new_seq;seq_bloc];
                end
            end
        end
        S(i,:) = new_seq;
    end
end
