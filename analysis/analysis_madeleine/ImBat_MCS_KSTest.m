% K-S Test to see if the join probabilities are sig diff.
function [h,p,h_row1,p_row1] = ImBat_MCS_KSTest(Joint_sim_Fnorm_vec,Joint_OG_Fnorm_vec,FnormM,OG_Fnorm)

    clear p h;
    [h,p] = kstest2(Joint_sim_Fnorm_vec,Joint_OG_Fnorm_vec);
    [h_row1,p_row1] = kstest2(FnormM(1,:),OG_Fnorm(1,:));

end

