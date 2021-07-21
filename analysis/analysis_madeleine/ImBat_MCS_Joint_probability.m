function [Joint_OG_Fnorm_vec,Joint_sim_Fnorm_vec] = ImBat_MCS_Joint_Probability(OG_Fnorm,FnormM,c_s_34,S)

% Calculate the joint probability of each flight-->flight pair

%% Real Data
% Make joint probability T by mulilpying each row by the probability of i
Joint_OG_Fnorm = [];
for i=1:size(OG_Fnorm,1)
    p_i = size(c_s_34(c_s_34==i),1)/size(c_s_34,1);
    Joint_OG_Fnorm(i,:) = OG_Fnorm(i,:)*p_i;
end
Joint_OG_Fnorm_vec = reshape(Joint_OG_Fnorm,1,size(Joint_OG_Fnorm,1)^2);

%% Simulated Data
% Make joint probability T by mulilpying each row by the probability of i

% Do the same for simulated data
% Get prior on each flight
for i=1:size(FnormM,1)
    for j=1:size(S,1)
        Srow = S(j,:);
        S_count(j) = size(Srow(Srow==i),2);
    end
    S_count_mean = mean(S_count);
    P_I(i,:) = S_count_mean/size(c_s_34,1);
end

Joint_sim_Fnorm = [];
for i=1:size(FnormM,1)
    p_i = size(c_s_34(c_s_34==i),1)/size(c_s_34,1);
    Joint_sim_Fnorm(i,:) = FnormM(i,:)*P_I(i);
end
Joint_sim_Fnorm_vec = reshape(Joint_sim_Fnorm,1,size(FnormM,1)^2);

