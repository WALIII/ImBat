
function ImBat_GLM_fit(FlightAlignedROI);
% estimate the tuning of individual neurons and the population using GLM




clear S2save F2save other
%[S2save,F2save,other] =  ImBat_Spikes(FlightAlignedROI);
disp('aligning data');
[S2save,F2save,other] =  ImBat_Spikes(FlightAlignedROI);

% [S2save,F2save] =  ImBat_Spikes(FlightAlignedROI_rando2);

% S2save(S2save==0) = NaN;

% calculate:

%%% Velocity
clear x1 x2 y1 y2 z1 z2 s spd
x1 = F2save(2:end,1);
x2 = F2save(1:end-1,1);
y1 = F2save(2:end,2);
y2 = F2save(1:end-1,2);
z1 = F2save(2:end,3);
z2 = F2save(1:end-1,3);

S=sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);
S(S>40) = 0;
spd = squeeze(S);

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
X_r = round(downsample(other.T_reward(1:end-1),ds_rate));
X_time = cat(2,round(downsample(other.T_forward(1:end-1),ds_rate)/ds_place),round(downsample(other.T_reverse(1:end-1),ds_rate)/ds_place));
X_theta = cat(2,round(downsample(theta,ds_rate)),round(downsample(theta2,ds_rate)),round(downsample(theta3,ds_rate)));

X_dist = round(downsample(other.T_dist(1:end-1),ds_rate));


X_xyz = cat(2,round(downsample(x1,ds_rate)/ds_place),round(downsample(y1,ds_rate)/ds_place),round(downsample(z1,ds_rate)/ds_place),X_theta);

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
mdl_time = fitglm(X_time,y,'Distribution','poisson');
mdl_v = fitglm(X_v,y,'Distribution','poisson');
mdl_r = fitglm(X_r,y,'Distribution','poisson');
mdl_theta = fitglm(X_theta,y,'Distribution','poisson');

mdl_all = fitglm(X,y,'Distribution','poisson');


p_vel(counter) = mdl_v.Rsquared.Adjusted;
p_xyz(counter) = mdl_xyz.Rsquared.Adjusted;
p_time(counter) = mdl_time.Rsquared.Adjusted;
p_r(counter) = mdl_r.Rsquared.Adjusted;
p_theta(counter) = mdl_theta.Rsquared.Adjusted;

p_all(counter) = mdl_all.Rsquared.Adjusted;

    counter = counter+1
    catch
    end
end

% remove sub zero scores
p_vel(p_vel<0) = 0;
p_all(p_all<0) = 0;
p_xyz(p_xyz<0) = 0;
p_time(p_time<0) = 0;
p_r(p_r<0) = 0;



% Plotting:
figure();
hold on;
histogram(p_vel,'BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs');
histogram(p_all,'BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs');
% histogram(p_vel,'BinWidth',0.01,'Normalization','probability');
% histogram(p_position,'BinWidth',0.01,'Normalization','probability');

figure();
boxplot([p_vel' p_time'  p_theta' p_xyz' p_all'],'Labels',{'vel','time', 'angle','pos', 'all'},'Notch','on');
ylabel(' R^2')
title('Spiking variance explained given predictors')






% figure();
% hold on;
% for i = 1:size(p_vel,2);
% x = [ 1 2];
% y = [p_vel(i) p_all(i)];;
% plot(x, y, '-k')
% hold on
% scatter(x(1), y(1), 50, 'b', 'filled')
% scatter(x(2), y(2), 50, 'r', 'filled')
%  
% end
%1-var(stats.resid)./var(y);