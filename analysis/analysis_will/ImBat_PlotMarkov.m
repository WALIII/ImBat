function ImBat_PlotMarkov(out_markov,FlightAlignedROI,ROI2use);


% plot output from:
% [out_markov] = ImBat_New_Markov(flightPaths);

% WAL3;
% d11.25.2020;
ROI_ON = 120

% get the clust number from aligned ROI..

clust_number = FlightAlignedROI.clust_number;
Ref_Clust = 5; % check this transition...




% colors:
col = hsv(4);


% plotting:
plt1 = 1;

% Find transitions in:
% 1-->A
% 2-->A
% 3-->A
% 4-->A
% 5-->A

% Find transitions out:
% A-->1
% A-->2
% A-->3
% A-->4
% A-->5



VA = out_markov.VA;
sort2use = out_markov.out_sort;

A = clust_number;
B = Ref_Clust;

% A-->1
a_1 = strfind(VA',[A 1]);
% A-->2
a_2 = strfind(VA',[A 2]); % add one to index into A...
% A-->3
a_3 = strfind(VA',[A 3]);
% A-->4
a_4 = strfind(VA',[A 4]);
% A-->5
a_5 = strfind(VA',[A 5]);
% A-->6
a_6 = strfind(VA',[A 6]);


% Compare to:

% A-->1
f1_a = strfind(VA',[1 A])+1;
% A-->2
f2_a = strfind(VA',[2 A])+1; % add one to index into A...
% A-->3
f3_a = strfind(VA',[3 A])+1;
% A-->4
f4_a = strfind(VA',[4 A])+1;
% A-->5
f5_a = strfind(VA',[5 A])+1;
% A-->6
f6_a = strfind(VA',[6 A])+1;

% idx = ~ismember(a_all,a_a); % find if a-->b transition
% a_x(idx==0) = []; % remove A-->B transitions
% 
% % X-->A
% idx2 = ~ismember(a_all,b_a); % find if b-->a transition
% x_a(idx2==0) = []; % remove B-->A transitions


% See if flights look different:
idX{1} = a_1;
idX{2} = a_2; 
idX{3} = a_3; 
idX{4} = a_4;
idX{5} = a_5;
idX{6} = a_6;


idX2{1} = f1_a;
idX2{2} = f2_a;
idX2{3} = f3_a;
idX2{4} = f4_a;
idX2{5} = f5_a;
idX2{6} = f6_a;



ROIMAT = zscore(squeeze(FlightAlignedROI.C_raw(ROI2use,:,:)));

titchr_1 = {'A-->1','A-->2','A-->3','A-->4','A-->5','A-->6'}

figure();
for i = 1:size(idX,2);
    for ii = 1: size(idX{i},2);
    a(ii) = find(FlightAlignedROI.cluster_idX == sort2use(idX{i}(ii)));
    ROImat(:,ii) = ROIMAT(:,a(ii));
    end

    try

   subplot(1,size(idX2,2),i)
aaa = sum(ROImat);
aaa2rmv = find(aaa ==0);
ROImat(:,aaa2rmv) = [];
clear aaa aaa2rmv
   imagesc(ROImat',[-2 4]);
       % Get axis handle
       
    ax = gca;
    % Set where ticks will be
    ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
    % Set TickLabels;
    ylabel('filghts')
    xlabel('time from takeoff');
    ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};
   ROImat_pre{i} = ROImat;

    catch
        ROImat_pre{i} = 0;
    end
                title(titchr_1{i});

    clear ROImat;
   clear a;
end


titchr_2 = {'1-->A','2-->A','3-->A','4-->A','5-->A','6-->A'}
figure();
for i = 1:size(idX2,2);
    for ii = 1: size(idX2{i},2);
    a(ii) = find(FlightAlignedROI.cluster_idX == sort2use(idX2{i}(ii)));
    ROImat(:,ii) = ROIMAT(:,a(ii));
    end
    try
   subplot(1,size(idX2,2),i)
   aaa = sum(ROImat);
aaa2rmv = find(aaa ==0);
ROImat(:,aaa2rmv) = [];
clear aaa aaa2rmv
   imagesc(ROImat',[-2 4]);
       % Get axis handle
    ax = gca;
    % Set where ticks will be
    ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
    % Set TickLabels;
    ylabel('filghts')
    xlabel('time from takeoff');
    ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};
   ROImat_post{i} = ROImat;

    catch
           ROImat_post{i} = 0;
    end
      title(titchr_2{i});

    clear ROImat;
   clear a;
end


figure();
   aaa = sum(ROIMAT);
aaa2rmv = find(aaa ==0);
ROIMAT(:,aaa2rmv) = [];
clear aaa aaa2rmv
imagesc(ROIMAT(:,:)',[-2 4]);
    % Get axis handle
    ax = gca;
    % Set where ticks will be
    ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
    % Set TickLabels;
    ylabel('filghts')
    xlabel('time from takeoff');
    ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};
    col = jet(size(idX2,2)+1);

figure(); 
subplot(1,2,1);
hold on;
title('pre');
for i = 1:size(idX2,2);
adata = ROImat_pre{i}';
cc = sum(adata,2);
ind2rmv = find(cc==0);
adata(ind2rmv,:) = [];
clear ind2rmv cc
if(size(adata,1)>4);
L = size(adata,2);
se = std(adata)/2;%sqrt(length(adata));
mn = median(adata);
%mn = mn-mean(mn(1:200));
mn = mn-mean(mn(650:750));
mn = smooth(mn,20)';
try
h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col(i,:),'DisplayName',titchr_1{i}); alpha(0.5);
plot(mn,'Color',col(i,:),'DisplayName','');
catch
end
legend('show'); %create/show legend

end
end


subplot(1,2,2);
hold on;
title('pre');
for i = 1:size(idX2,2);
adata = ROImat_post{i}';
cc = sum(adata,2);
ind2rmv = find(cc==0);
adata(ind2rmv,:) = [];
clear ind2rmv cc
if(size(adata,1)>4);
L = size(adata,2);
se = std(adata)/2;%sqrt(length(adata));
mn = median(adata);
%mn = mn-mean(mn(1:200));
mn = mn-mean(mn(650:750));
mn = smooth(mn,20)';
h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col(i,:),'DisplayName',titchr_2{i}); alpha(0.5);
plot(mn,'Color',col(i,:),'DisplayName','');
% legend(titchr_2{i},'');
end
legend('show'); %create/show legend

end


% 
% % STATs
% 
% adata = ROImat_post{2}';
% adata2 = ROImat_post{4}';
% x1 = corr(cat(1,adata,adata2)');
% 
% x1_a = corr(adata');
% x1_a = triu(x1_a); x1_a = x1_a(:); x1_a(x1_a==1) = []; x1_a(x1_a==0) = []; 
% x2_a = corr(adata2');
% x2_a = triu(x2_a); x2_a = x2_a(:); x2_a(x2_a==1) = []; x2_a(x2_a==0) = []; 
% 
% x1_nonOverlap = cat(1,x2_a,x1_a);
% x1_overlap = x1(size(adata,1):end,1:size(adata,1)); x1_overlap = x1_overlap(:);
% 
% [pval_combined_data,~] = ranksum(((x1_overlap)), ((x1_nonOverlap)))
% 
% 
% figure(); hold on;
% histogram(x1_nonOverlap,'Normalization','probability');
% histogram(x1_overlap,'Normalization','probability');
% 
% [pval_combined_data,~] = ranksum(((dist1)), ((dist2)))
% end

