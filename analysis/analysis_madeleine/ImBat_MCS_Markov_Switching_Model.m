% Markov Switching Model

% Fit an arima to the beginning of each session
[ss,rr] = sort(flightPaths34.flight_starts_idx);
day_idx = flightPaths34.day(rr);

Bs = {};
beginnings = [];
middles_1 = []; middles_2 = [];
for j=1:length(unique(day_idx))
    Bs{j} = find(day_idx == j)
end
for i=1:length(Bs)
    temp = c_s_34(Bs{i});
    beginnings = [beginnings,temp(1:10)];
    if length(temp) > 25
        middles_1 = [middles_1,temp(15:25)];
    end
    if length(temp) > 50
        middles_2 = [middles_2,temp(end-35:end-15)];
    end
end
Beginnings_training = reshape(beginnings,size(beginnings,1)*size(beginnings,2),1);
Middles_training = reshape(middles_2,size(middles_2,1)*size(middles_2,2),1);

% Template Model (Unclustered)
Mdl_1 = arima(2,0,1);
EstMdl_1 = estimate(Mdl,Beginnings_training);
EstMdl_1.Description = "Unclustered";

% Template Model (Clustered)
Mdl_2 = arima(2,0,1);
EstMdl_2 = estimate(Mdl,Middles_training);
EstMdl_1.Description = "Clustered";


% Make Markov Switching Model
P_clustered_to_unclustered = c_s_Tnorm_34;
mc = dtmc(P_clustered_to_unclustered,'StateNames',["Unclustered" "Clustered"])

mdl = [EstMdl_1; EstMdl_2];

% Partially specified markov switch model
Mdl = msVAR(mc,mdl);

% Model to estimate
mc0 = dtmc(P_clustered_to_unclustered,'StateNames',Mdl.StateNames);
mdl01 = arima(2,0,0);
mdl02 = arima(2,0,0);
Mdl_est = msVAR(mc0,[mdl01; mdl02]);

%Estimate 
EstMdl = estimate(Mdl_est,Mdl,c_s_34,'IterationPlot',true);