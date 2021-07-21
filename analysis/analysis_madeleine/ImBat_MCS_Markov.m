function Markov_Output = ImBat_MCS_Markov(flightPaths)


% MCS 3/31/21

% Path to data: /Users/madeleinesnyder/Desktop/Processed/CellReg_files/ROI_Data

% Script for building a Markov Model of the flight behavior.
% First-order transition statistics only in this version
% Based on https://www.mathworks.com/help/stats/hidden-markov-models-hmm.html#f8288

num_clust = size(flightPaths.clusterIndex,2);

% Create vector of all flight cluster IDs 

true_times = flightPaths.AllFlightsMasterTime(flightPaths.flight_starts_idx(:));
[tts,tts_idx]  = sort(true_times(:));
c_s = flightPaths.id(tts_idx);

% Get the state transition probabilities

nxt_clst = []; % each column is the vector of which flight that column index goes to
prv_clst = [];
for i=1:num_clust
    i_idx = find(c_s(:) == i);
    for j=1:size(i_idx,1)
        if i_idx(j) == size(c_s,1)
            prv_clst(i,j) = c_s(i_idx(j)-1);
        elseif i_idx(j) == 1
            nxt_clst(i,j) = c_s(i_idx(j)+1);
        else
            nxt_clst(i,j) = c_s(i_idx(j)+1);
            prv_clst(i,j) = c_s(i_idx(j)-1);
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
Fnorm = []; FPnorm = [];
for i=1:num_clust
    row = F(i,:);
    rowP = FP(i,:);
    Fnorm(i,:) = row./sum(row);
    FPnorm(i,:) = rowP./sum(rowP);
end

Fnorm_first5 = Fnorm(1:5,1:5);
Fnorm_noone = Fnorm(2:10,2:10);
Fnorm_first10 = Fnorm(1:10,1:10);
FPnorm_first10 = FPnorm(1:10,1:10);

figure(); title("Flight Cluster Transition Probability Matrix"); h = heatmap(Fnorm, 'Colormap', copper);
figure(); title("Flight Cluster Transition Probability Matrix"); h = heatmap(Fnorm_first10, 'Colormap', copper);
figure(); title("Flight Cluster Transition Probability Matrix History"); h = heatmap(FPnorm_first10, 'Colormap', copper);
figure(); title("Flight Cluster Transition Probability Matrix History"); h = heatmap(FPnorm, 'Colormap', copper);

% DTMC
mc10 = dtmc(Fnorm_first10); 
figure(); 
h = graphplot(mc10);
c = h.EdgeColor;
h.EdgeColor = 'k';
mcp = nonzeros(reshape(mc10.P,[size(mc10.P,1)*size(mc10.P,1),1])).*5;
h.LineWidth = mcp;
h.MarkerSize = [size(flightPaths.id(flightPaths.id==1),1)/10,
                size(flightPaths.id(flightPaths.id==2),1)/10,
                size(flightPaths.id(flightPaths.id==3),1)/10,
                size(flightPaths.id(flightPaths.id==4),1)/10,
                size(flightPaths.id(flightPaths.id==5),1)/10,
                size(flightPaths.id(flightPaths.id==6),1)/10,
                size(flightPaths.id(flightPaths.id==7),1)/10,
                size(flightPaths.id(flightPaths.id==8),1)/10,
                size(flightPaths.id(flightPaths.id==9),1)/10,
                size(flightPaths.id(flightPaths.id==10),1)/10]

% DTMC
mc = dtmc(Fnorm_first5); 
figure(); 
h = graphplot(mc);
c = h.EdgeColor;
h.EdgeColor = 'k';
mcp = nonzeros(reshape(mc.P,[size(mc.P,1)*size(mc.P,1),1])).*5;
h.LineWidth = mcp;
h.MarkerSize = [size(flightPaths.id(flightPaths.id==1),1)/10,
                size(flightPaths.id(flightPaths.id==2),1)/10,
                size(flightPaths.id(flightPaths.id==3),1)/10,
                size(flightPaths.id(flightPaths.id==4),1)/10,
                size(flightPaths.id(flightPaths.id==5),1)/10]
            
% DTMC exclude 1
mc10 = dtmc(Fnorm_noone); 
figure(); 
h = graphplot(mc10);
c = h.EdgeColor;
h.EdgeColor = 'k';
mcp = nonzeros(reshape(mc10.P,[size(mc10.P,1)*size(mc10.P,1),1])).*5;
h.LineWidth = mcp;
h.MarkerSize = [size(flightPaths.id(flightPaths.id==2),1)/10,
                size(flightPaths.id(flightPaths.id==3),1)/10,
                size(flightPaths.id(flightPaths.id==4),1)/10,
                size(flightPaths.id(flightPaths.id==5),1)/10,
                size(flightPaths.id(flightPaths.id==6),1)/10,
                size(flightPaths.id(flightPaths.id==7),1)/10,
                size(flightPaths.id(flightPaths.id==8),1)/10,
                size(flightPaths.id(flightPaths.id==9),1)/10,
                size(flightPaths.id(flightPaths.id==10),1)/10]

            
% DTMC ALL
mc10 = dtmc(Fnorm); 
figure(); 
h = graphplot(mc10);
c = h.EdgeColor;
h.EdgeColor = 'k';
mcp = nonzeros(reshape(mc10.P,[size(mc10.P,1)*size(mc10.P,1),1])).*5;
h.LineWidth = mcp;
h.MarkerSize = [size(flightPaths.id(flightPaths.id==1),1)/10,
                size(flightPaths.id(flightPaths.id==2),1)/10,
                size(flightPaths.id(flightPaths.id==3),1)/10,
                size(flightPaths.id(flightPaths.id==4),1)/10,
                size(flightPaths.id(flightPaths.id==5),1)/10,
                size(flightPaths.id(flightPaths.id==6),1)/10,
                size(flightPaths.id(flightPaths.id==7),1)/10,
                size(flightPaths.id(flightPaths.id==8),1)/10,
                size(flightPaths.id(flightPaths.id==9),1)/10,
                size(flightPaths.id(flightPaths.id==10),1)/10,
                size(flightPaths.id(flightPaths.id==11),1)/10,
                size(flightPaths.id(flightPaths.id==12),1)/10,
                size(flightPaths.id(flightPaths.id==13),1)/10,
                size(flightPaths.id(flightPaths.id==14),1)/10,
                size(flightPaths.id(flightPaths.id==15),1)/10,
                size(flightPaths.id(flightPaths.id==16),1)/10,
                size(flightPaths.id(flightPaths.id==17),1)/10,
                size(flightPaths.id(flightPaths.id==18),1)/10,
                size(flightPaths.id(flightPaths.id==19),1)/10,
                size(flightPaths.id(flightPaths.id==20),1)/10,
                size(flightPaths.id(flightPaths.id==21),1)/10,
                size(flightPaths.id(flightPaths.id==22),1)/10,
                size(flightPaths.id(flightPaths.id==23),1)/10,
                size(flightPaths.id(flightPaths.id==24),1)/10,
                size(flightPaths.id(flightPaths.id==25),1)/10,
                size(flightPaths.id(flightPaths.id==26),1)/10,
                size(flightPaths.id(flightPaths.id==27),1)/10,
                size(flightPaths.id(flightPaths.id==28),1)/10,
                size(flightPaths.id(flightPaths.id==29),1)/10,
                size(flightPaths.id(flightPaths.id==30),1)/10,
                size(flightPaths.id(flightPaths.id==31),1)/10,
                size(flightPaths.id(flightPaths.id==32),1)/10,
                size(flightPaths.id(flightPaths.id==33),1)/10,
                size(flightPaths.id(flightPaths.id==34),1)/10,
                size(flightPaths.id(flightPaths.id==35),1)/10]     
            

figure(); 
subplot(10,1,1); bar(Fnorm_first10(1,:)); ylim([0,1]);
subplot(10,1,2); bar(Fnorm_first10(2,:)); ylim([0,1]);
subplot(10,1,3); bar(Fnorm_first10(3,:)); ylim([0,1]);
subplot(10,1,4); bar(Fnorm_first10(4,:)); ylim([0,1]);
subplot(10,1,5); bar(Fnorm_first10(5,:)); ylim([0,1]);
subplot(10,1,6); bar(Fnorm_first10(6,:)); ylim([0,1]);
subplot(10,1,7); bar(Fnorm_first10(7,:)); ylim([0,1]);
subplot(10,1,8); bar(Fnorm_first10(8,:)); ylim([0,1]);
subplot(10,1,9); bar(Fnorm_first10(9,:)); ylim([0,1]);
subplot(10,1,10); bar(Fnorm_first10(10,:)); ylim([0,1]);


figure(); 
subplot(10,1,1); bar(Fnorm(1,1:10)); ylim([0,1]);
subplot(10,1,2); bar(Fnorm(2,:)); ylim([0,1]);
subplot(10,1,3); bar(Fnorm(3,:)); ylim([0,1]);
subplot(10,1,4); bar(Fnorm(4,:)); ylim([0,1]);
subplot(10,1,5); bar(Fnorm(5,:)); ylim([0,1]);
subplot(10,1,6); bar(Fnorm(6,:)); ylim([0,1]);
subplot(10,1,7); bar(Fnorm(7,:)); ylim([0,1]);
subplot(10,1,8); bar(Fnorm(8,:)); ylim([0,1]);
subplot(10,1,9); bar(Fnorm(9,:)); ylim([0,1]);
subplot(10,1,10); bar(Fnorm(10,:)); ylim([0,1]);

figure(); 
subplot(10,1,1); bar(Fnorm(11,:)); ylim([0,1]);
subplot(10,1,2); bar(Fnorm(12,:)); ylim([0,1]);
subplot(10,1,3); bar(Fnorm(13,:)); ylim([0,1]);
subplot(10,1,4); bar(Fnorm(14,:)); ylim([0,1]);
subplot(10,1,5); bar(Fnorm(15,:)); ylim([0,1]);
subplot(10,1,6); bar(Fnorm(16,:)); ylim([0,1]);
subplot(10,1,7); bar(Fnorm(17,:)); ylim([0,1]);
subplot(10,1,8); bar(Fnorm(18,:)); ylim([0,1]);
subplot(10,1,9); bar(Fnorm(19,:)); ylim([0,1]);
subplot(10,1,10); bar(Fnorm(20,:)); ylim([0,1]);


% Look at transitions between stereotyped and non-stereotyped
c_s_stereo = [];
for i=1:size(c_s,1)
    if c_s(i) == 1
        c_s_stereo(i) = 2;
    elseif c_s(i) < 10
        c_s_stereo(i) = 1;
    end
end

% Calculate transition proability between 1 (stereotyped and 2
% (understereotyped)
c_s_T = zeros(2,2); 
for i=1:size(c_s,1)
    if i == size(c_s,1)
        disp('')
        break
    end
    if c_s_stereo(i) == 1 & c_s_stereo(i+1) == 1
        c_s_T(1,1) = c_s_T(1,1)+1;
    elseif c_s_stereo(i) == 1 & c_s_stereo(i+1) == 2
        c_s_T(1,2) = c_s_T(1,2)+1;
    elseif c_s_stereo(i) == 2 & c_s_stereo(i+1) == 1
        c_s_T(2,1) = c_s_T(2,1)+1;
    elseif c_s_stereo(i) == 2 & c_s_stereo(i+1) == 2
        c_s_T(2,2) = c_s_T(2,2)+1;
    end
end

c_s_Tnorm_first10(1,:) = c_s_T(1,:)./sum( c_s_T(1,:));
c_s_Tnorm_first10(2,:) = c_s_T(2,:)./sum( c_s_T(2,:));

% Probability of transitioning from 
% (stereotyped to stereotyped)
%                               (unstereotyped to unstereotyped)
% Including the first 10 clusters in the stereotyped group

%    0.8898    0.1102
%    0.3171    0.6829


% DTMC
mc = dtmc(Fnorm_first5); 
figure(); 
h = graphplot(mc);
c = h.EdgeColor;
h.EdgeColor = 'k';
mcp = nonzeros(reshape(mc.P,[size(mc.P,1)*size(mc.P,1),1])).*5;
h.LineWidth = mcp;
h.MarkerSize = [size(flightPaths.id(flightPaths.id==1),1)/10,
                size(flightPaths.id(flightPaths.id==2),1)/10,
                size(flightPaths.id(flightPaths.id==3),1)/10,
                size(flightPaths.id(flightPaths.id==4),1)/10,
                size(flightPaths.id(flightPaths.id==5),1)/10];




end
   