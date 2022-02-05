
function output = ImBat_GLM_fit(FlightAlignedROI);
% estimate the tuning of individual neurons and the population using GLM




clear S2save F2save other
%[S2save,F2save,other] =  ImBat_Spikes(FlightAlignedROI);
disp('aligning data');
[S2save,F2save,other] =  ImBat_Spikes(FlightAlignedROI);

% [S2save,F2save] =  ImBat_Spikes(FlightAlignedROI_rando2);

% S2save(S2save==0) = NaN;

% calculate:

%%% Velocity
% clear x1 x2 y1 y2 z1 z2 s spd
x1 = F2save(2:end,1);
x2 = F2save(1:end-1,1);
y1 = F2save(2:end,2);
y2 = F2save(1:end-1,2);
z1 = F2save(2:end,3);
z2 = F2save(1:end-1,3);
% 
% S=sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);
% S(S>40) = 0;
spd = other.T_spd;

theta = zeros(1, size(x1,1));
theta2 = zeros(1, size(x1,1));
theta3 = zeros(1, size(x1,1));

%%% angle (azmuthal)
for i = 1:size(x1,1)-60
    theta(i) = atan2(y1(i+60)-y1(i),x1(i+60)-x1(i));
    theta2(i) = atan2(y1(i+60)-y1(i),z1(i+60)-z1(i));
    theta3(i) = atan2(z1(i+60)-z1(i),x1(i+60)-x1(i));
end
theta = theta';
theta2 = theta2';
theta3 = theta3';


% Time from takeoff
% Distance from Landing
ds_rate = 20;
ds_place = 10;
X_v = round(downsample(spd,ds_rate));
X_dist = round(downsample(other.T_dist(1:end-1),ds_rate));

X_r = round(downsample(other.T_reward(1:end-1),ds_rate));
X_time = cat(2,round(downsample(other.T_forward(1:end-1),ds_rate)/ds_place),round(downsample(other.T_reverse(1:end-1),ds_rate)/ds_place));
X_theta = cat(2,round(downsample(theta,ds_rate)),round(downsample(theta2,ds_rate)),round(downsample(theta3,ds_rate)));
X_xyz = cat(2,round(downsample(x1,ds_rate)/ds_place),round(downsample(y1,ds_rate)/ds_place),round(downsample(z1,ds_rate)/ds_place));

% Combinations:
X_xyz_theta = cat(2,round(downsample(x1,ds_rate)/ds_place),round(downsample(y1,ds_rate)/ds_place),round(downsample(z1,ds_rate)/ds_place),X_theta);
X_time_xyz = cat(2,X_time,X_xyz);

% Tripple Combinations:
X_time_xyz_vel = cat(2,X_time,X_xyz,X_v); % no theta
X_time_xyz_theta = cat(2,X_time,X_xyz,X_theta); % no velocity
X_time_theta_vel = cat(2,X_time,X_theta,X_v); % no position
X_xyz_theta_vel = cat(2,X_xyz,X_theta,X_v); % no time


X = cat(2,round(downsample(x1,ds_rate)/ds_place),round(downsample(y1,ds_rate)/ds_place),round(downsample(z1,ds_rate)/ds_place),round(downsample(spd,ds_rate)),X_time,X_theta,X_dist);

disp('Fitting models to data');
clear p_position p_vel
counter = 1;
warning off
for i = 1:100;
    try
        y = S2save(1:end-1,i);
        y = smooth(y,25);
        y = ((downsample(y,ds_rate)));
        y(y==0) = NaN;
        mdl_xyz = fitglm(X_xyz,y,'Distribution','poisson');
        mdl_xyz_theta = fitglm(X_xyz_theta,y,'Distribution','poisson');

        mdl_time = fitglm(X_time,y,'Distribution','poisson');
        mdl_time_xyz = fitglm(X_time_xyz,y,'Distribution','poisson');
        mdl_v = fitglm(X_v,y,'Distribution','poisson');
        mdl_r = fitglm(X_r,y,'Distribution','poisson');
        mdl_theta = fitglm(X_theta,y,'Distribution','poisson');
        mdl_dist = fitglm(X_dist,y,'Distribution','poisson');
        
        mdl_all = fitglm(X,y,'Distribution','poisson');
        
        % combinations
        mdl_time_xyz_vel = fitglm(X_time_xyz_vel,y,'Distribution','poisson');
        mdl_time_xyz_theta = fitglm(X_time_xyz_theta,y,'Distribution','poisson');% no velocity
        mdl_time_theta_vel = fitglm(X_time_theta_vel,y,'Distribution','poisson');% no position
        mdl_xyz_theta_vel = fitglm(X_xyz_theta_vel,y,'Distribution','poisson');% no time
       
    
        p_vel(counter) = mdl_v.Rsquared.Adjusted;
        p_xyz(counter) = mdl_xyz.Rsquared.Adjusted;
        p_time(counter) = mdl_time.Rsquared.Adjusted;
        p_time_xyz(counter) = mdl_time_xyz.Rsquared.Adjusted;
        p_r(counter) = mdl_r.Rsquared.Adjusted;
        p_theta(counter) = mdl_theta.Rsquared.Adjusted;
        p_dist(counter) = mdl_dist.Rsquared.Adjusted;
        p_xyz_theta(counter) = mdl_xyz_theta.Rsquared.Adjusted;
        p_all(counter) = mdl_all.Rsquared.Adjusted;
        
        % Triple combinations
        
        p_time_xyz_vel(counter) = mdl_time_xyz_vel.Rsquared.Adjusted; % no theta
        p_time_xyz_theta(counter) = mdl_time_xyz_theta.Rsquared.Adjusted; % no velocity
        p_time_theta_vel(counter) = mdl_time_theta_vel.Rsquared.Adjusted; % no position
        p_xyz_theta_vel(counter) = mdl_xyz_theta_vel.Rsquared.Adjusted; % no time
       
        
        Unique_theta(counter) = p_all(counter)-p_time_xyz_vel(counter);
        Unique_vel(counter) = p_all(counter)-p_time_xyz_theta(counter);
        Unique_xyz(counter) = p_all(counter)-p_time_theta_vel(counter);
        Unique_time(counter) = p_all(counter)-p_xyz_theta_vel(counter);

        % now get pvalues:
       pVal(1,counter)   = mdl_v.devianceTest.pValue(2);
       pVal(2,counter)   = mdl_xyz.devianceTest.pValue(2); 
       pVal(3,counter)   = mdl_time.devianceTest.pValue(2); 
       pVal(4,counter)   = mdl_r.devianceTest.pValue(2); 
       pVal(5,counter)   = mdl_theta.devianceTest.pValue(2); 
       pVal(6,counter)   = mdl_dist.devianceTest.pValue(2); 

       
       pVal(7,counter)   = mdl_time_xyz.devianceTest.pValue(2); 
       pVal(8,counter)   = mdl_xyz_theta.devianceTest.pValue(2); 

       pVal(9,counter)   = mdl_time_xyz_vel.devianceTest.pValue(2); 
       pVal(10,counter)  = mdl_all.devianceTest.pValue(2); 

        counter = counter+1;
    catch
        disp('huh?');
    end
    
    
end

% remove sub zero scores
p_vel(p_vel<0) = 0;
% p_all(p_all<0) = 0;
p_xyz(p_xyz<0) = 0;
p_xyz_theta(p_xyz_theta<0) = 0;

p_time(p_time<0) = 0;
p_time_xyz(p_time_xyz<0) = 0;
p_r(p_r<0) = 0;
p_dist(p_dist<0) = 0;

% remove zero for combined:
%Unique_time(Unique_time<0) = 0;

% Plotting:
figure();
hold on;
histogram(p_vel,'BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs');
histogram(p_all,'BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs');
% histogram(p_vel,'BinWidth',0.01,'Normalization','probability');
% histogram(p_position,'BinWidth',0.01,'Normalization','probability');


figure();
boxplot([abs(Unique_time') abs(Unique_theta') abs(Unique_vel') abs(Unique_xyz') abs(p_all')],'Labels',{'unique_time', 'unique_theta', 'unique_vel', 'unique_xyz', 'Combined'},'Notch','on');
ylabel(' R^2')
title('Spiking variance explained given predictors')

[pval_combined_data,~] = ranksum(abs(Unique_xyz'),abs(p_all'));

[pval_xyz_vs_time,~] = ranksum(abs(Unique_xyz'),abs(Unique_time'))
[pval_xyz_vs_theta,~] = ranksum(abs(Unique_xyz'),abs(Unique_theta'))
[pval_xyz_vs_velocity,~] = ranksum(abs(Unique_xyz'),abs(Unique_vel'))


figure();
boxplot([p_vel' p_time'  p_dist' p_theta' p_xyz' p_xyz_theta' p_time_xyz' p_all'],'Labels',{'vel','time', 'distance','theta','pos','theta+pos', 'time+pos','all'},'Notch','on');
ylabel(' R^2')
title('Spiking variance explained given predictors')

% venn diagram
scfact = 2; 
area1 = mean(p_time)*scfact; 
area2 = mean(p_xyz)*scfact;
area3 = mean(p_vel)*scfact;

overlap1 = mean(p_time_xyz_vel) - ( (mean(p_time_xyz_vel)-mean(p_xyz)) + (mean(p_time_xyz_vel)-mean(p_vel)))/2;
overlap2 = mean(p_time_xyz_vel) - ( (mean(p_time_xyz_vel)-mean(p_time)) + (mean(p_time_xyz_vel)-mean(p_vel)))/2;
overlap3 = mean(p_time_xyz_vel) - ( (mean(p_time_xyz_vel)-mean(p_xyz)) + (mean(p_time_xyz_vel)-mean(p_time)))/2;
overlap4 = mean(p_time_xyz_vel) - ( (mean(p_time_xyz_vel)-mean(p_xyz)) + (mean(p_time_xyz_vel)-mean(p_time)) + (mean(p_time_xyz_vel)-mean(p_vel)))/2;

figure(); venn([area1 area2 area3], [overlap1 overlap2 overlap3 overlap4 ]);


% venn diagram
area1 = mean(p_time); 
area2 = mean(p_xyz);

overlap1 = mean(p_xyz)-(mean(p_time_xyz)-mean(p_time));
overlap2 = mean(p_time)-(mean(p_time_xyz)-mean(p_xyz));


figure(); venn([area1 area2], [overlap2 overlap2]);

% EXPORT DATA
output.pVal = pVal;

for i = 1:10;
    p2plot(i) = size(find(pVal(i,:)<0.001),2)/size(pVal,2)*100;
end

figure(); barh(p2plot);
somenames={'vel'; 'pos'; 'time'; 'reward'; 'head'; 'time +'; 'time+pos'; 'head+pos'; 'time+poss+vel'; 'Combined' };
set(gca,'yticklabel',somenames)
xlabel('Percent of cells significantly modulated')


% plot exclusivity
A = find(pVal(2,:)>0.001 & pVal(3,:)<0.001);
output.Time_exclusive = size(A,2)/size(pVal,2)*100
B = find(pVal(2,:)<0.001 & pVal(3,:)>0.001);
output.Place_exclusive = size(B,2)/size(pVal,2)*100


% find the 'nots'
