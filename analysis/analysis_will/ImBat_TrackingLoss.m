function ImBat_TrackingLoss(out_dat)

% plot output of across day ROI metrics
% Calculate tracking loss:

% NOTE manually loaded in all the ROI matrics:
% out_dat{1} = ImBat_analysis_12_01_2020(cell_registered_struct,ROI_Data)


% calculate: out of all cells, what proportion of x do we hold for n days?
% take the longest tracked interval for all cells, then compute:
% 1. total-(tracked interval)./total ( % cells lost, relative to total)
% 2. cummulative loss: ( survival fraction);


figure(); 
hold on;
for i = 1:size(out_dat,2);
    plot(out_dat{i}.tracked_ROIs);
end
title('ROIs/day');
ylabel('ROIs detected');
xlabel('days');


% get CDF of tracking stability:

hold on;
clear ff f_temp
for i = 1:size(out_dat,2);
    if i ==1;
        ff = sum(out_dat{i}.tracking_stability,2);
    else
        f_temp = sum(out_dat{i}.tracking_stability,2);
    ff = cat(1,ff,f_temp);
    clear f_temp
    end
end
[aa ab]  = sort(ff);
figure(); plot(aa, 1:length(aa));

for i = 1:12;
    a_s(i) = size(find(aa == i),1);
end

a_s(1) = [];

figure(); bar(a_s/size(aa,1)*100);
size(find(aa > 3),1)/size(aa,1)


% use time, instead of sessions:
% get CDF of tracking stability:

hold on;
clear ff f_temp
fS = zeros(1,40);
for i = 1:size(out_dat,2);
    if i ==1;
        % first, remove cells that are not tracked:
        ffdat = out_dat{i}.tracking_stability;
        ffdat(sum(ffdat,2)==1,:) = [];

        fS(out_dat{i}.days') = fS(out_dat{i}.days')+nansum(ffdat);
        ffdat = ffdat'.*out_dat{i}.days';
        
%         ffdat(ffdat==0) = NaN;
        ff = (max(ffdat)-min(ffdat))';
        ff2 = (max(sort(ffdat')')-min(sort(ffdat')'))';
        
        % get normalized:
        [aa ab]  = sort(ff);
for ii = 1:max(ff);
    a_s1(ii) = size(find(aa == ii),1);
end
A2save(i,:) = a_s1;
A2save_percent(i,:) = (sum(a_s1)-a_s1)/sum(a_s1);
SurvivalFraction(i,:) = -cumsum(-((sum(a_s1)-a_s1)/sum(a_s1))+1)+1;

clear aa a_s1


       

    else
        clear ffdat
           ffdat = out_dat{i}.tracking_stability;
        ffdat(sum(ffdat,2)==1,:) = [];

        fS(out_dat{i}.days') = fS(out_dat{i}.days')+nansum(ffdat);
        ffdat = ffdat'.*out_dat{i}.days';
        ffdat(ffdat==0) = NaN;
        f_temp = (max(ffdat)-min(ffdat))';
        ftemp = (max(sort(ffdat')')-min(sort(ffdat')'))';

    ff = cat(1,ff,f_temp);
    ff2 = cat(1,ff2, ftemp);
    
            % get normalized:
[aa ab]  = sort(f_temp);
for ii = 1:max(f_temp);
    a_s1(ii) = size(find(aa == ii),1);
end
% A2save(ii,:) = a_s1;
htemp = (sum(a_s1)-a_s1)/sum(a_s1);
htemp2 = a_s1;
try
A2save_percent(i,:) = htemp(1:length(A2save_percent));
A2save(i,:) = htemp2(1:length(A2save_percent));

catch
    A2save_percent(i,(1:length(htemp))) = htemp;
    A2save(i,(1:length(htemp2))) = htemp2;
end


clear aa a_s1 htemp htemp2


    clear f_temp ffdat ftemp
    end
end
clear aa ab a_s
[aa ab]  = sort(ff);
[aa2 ab2]  = sort(ff2);
figure(); 
hold on;

plot(aa, 1:length(aa));
plot(aa2, 1:length(aa2));

title('Total length of days cells are tracked');

for i = 1:max(ff);
    a_s(i) = size(find(aa == i),1);
    a_s2(i) = size(find(aa2 == i),1);

    a_s_norm(i) = (size(find(aa == i),1)./fS(:,i+1));
end
% a_s(1) = [];


figure(); bar(a_s/a_s2);
size(find(aa > 3),1)/size(aa,1)
title('Total length of days cells are tracked');
ylabel('# ROIs')
xlabel('Days Tracked');

figure(); 
A2save_percent(A2save_percent==0) = NaN;

figure()
A2save_percent(A2save_percent==1) = NaN;
plot(nanmean(A2save_percent),'o')


% ========== stats ( regression) ===========%

clear y x1 X
y = (nanmean(A2save_percent)');

% add one to the front
y = cat(1,1,y);
X1 = (1:length(y))';
x1 = ones(size(X1,1),1);




X = [x1 X1];    % Includes column of ones

[~,~,~,~,stats] = regress(y,X);
stats
X = [X1];
mdl = fitlm(X,y);
figure();
plotAdded(mdl)
y2 = y;
X2 = X1;
y2(isnan(y)) = [];
X2(isnan(y)) = [];


coefs = polyfit(X2, y2, 1)

figure()
plot(X2,y2,'o')

figure(); plot(cumsum(sum(A2save))); title('CDF of all cells'); xlabel('days'); ylabel('Unique ROIs')

