function [out_markov] = ImBat_Markov(ROI_Data,day2use,nclust)

% Segregate Flights:
%day2use = 12;
Sim = 0;

[out] =  ImBat_SegTrajectories(ROI_Data{1,day2use}.Alignment.out.flights,ROI_Data{1,day2use}.Alignment.out.Location_time,'nclusters',nclust);
close all

% Plot flights

figure(); 
hold on;
plot3(ROI_Data{1,day2use}.Alignment.out.flights(:,1),ROI_Data{1,day2use}.Alignment.out.flights(:,2),ROI_Data{1,day2use}.Alignment.out.flights(:,3),'k');% plot the flight trajectory in space

% plot all segregated trajectories
hold on;
for i = 1:size(out.flight_starts_indx,2)
plot3(ROI_Data{1,day2use}.Alignment.out.flights(out.flight_starts_indx(i):out.flight_ends_indx(i),1),ROI_Data{1,day2use}.Alignment.out.flights(out.flight_starts_indx(i):out.flight_ends_indx(i),2),ROI_Data{1,day2use}.Alignment.out.flights(out.flight_starts_indx(i):out.flight_ends_indx(i),3),'r');% plot the flight trajectory in space
end


% Plot specific flights
figure();
hold on;
colorC = hsv(size(out.ClusterIndex,2));
for i = 1:size(out.ClusterIndex,2)
    for ii = 1: size(out.ClusterIndex{i},2)
        idX = out.ClusterIndex{i}(ii);
        plot3(ROI_Data{1,day2use}.Alignment.out.flights(out.flight_starts_indx(idX):out.flight_ends_indx(idX),1),ROI_Data{1,day2use}.Alignment.out.flights(out.flight_starts_indx(idX):out.flight_ends_indx(idX),2),ROI_Data{1,day2use}.Alignment.out.flights(out.flight_starts_indx(idX):out.flight_ends_indx(idX),3),'Color',colorC(i,:));% plot the flight trajectory in space
    end
end




% % plot all segregated trajectories
% hold on;
% for i = 1:size(out.flight_starts_indx,2)
% plot3(ROI_Data{1,day2use}.Alignment.out.flights(out.flight_starts_indx(i):out.flight_ends_indx(i),1),ROI_Data{1,day2use}.Alignment.out.flights(out.flight_starts_indx(i):out.flight_ends_indx(i),2),ROI_Data{1,day2use}.Alignment.out.flights(out.flight_starts_indx(i):out.flight_ends_indx(i),3),'r');% plot the flight trajectory in space
% end




% Plot specific flights
figure();
hold on;
colorC = hsv(size(out.ClusterIndex,2));
for i = 1:10%size(out.ClusterIndex,2)
    subplot(5,2,i);
    hold on;
    for ii = 1: size(out.ClusterIndex{i},2)
        idX = out.ClusterIndex{i}(ii);
        plot3(ROI_Data{1,day2use}.Alignment.out.flights(out.flight_starts_indx(idX):out.flight_ends_indx(idX),1),ROI_Data{1,day2use}.Alignment.out.flights(out.flight_starts_indx(idX):out.flight_ends_indx(idX),2),ROI_Data{1,day2use}.Alignment.out.flights(out.flight_starts_indx(idX):out.flight_ends_indx(idX),3),'Color',colorC(i,:));% plot the flight trajectory in space
    end
    title(['Flight ', num2str(i)]);
end




% build transition matrix:
Seq = zeros(1,size(out.flight_ends_indx,2));
for i = 1:size(out.ClusterIndex,2)
    Seq(out.ClusterIndex{i}) = i;
end


% build sequence
Seq2 = (1:size(out.flight_starts_indx,2));
VA(:,1) = Seq;
% VA(:,2) = Seq2;
 [u,~,n] = unique(VA,'rows');
 N = length(u);
 F = accumarray([n(1:end-1),n(2:end)],1,[N,N]);
 T = bsxfun(@rdivide,F,sum(F,2));

 
 % Make simple markov model
 mc = dtmc(T);
 figure;
graphplot(mc,'ColorEdges',true);

% Suppose that the initial state distribution is uniform. Compute the evolution of the distribution for 20 time steps.
numSteps = 20;
X = redistribute(mc,numSteps);
figure;
distplot(mc,X);

% Is reduceable 
tf = isreducible(mc)


% Simulate:
if Sim ==1;
 % Suppose that the initial state distribution is uniform. Compute the evolution of the distribution for 20 time steps.
numSteps = 20;
X = redistribute(mc,numSteps);
figure;
distplot(mc,X);
% Walk through simulation
numSteps = 50;
X = simulate(mc,numSteps);
figure(); simplot(mc,X,'FrameRate',0.5,'Type','graph');

% x0 = [0 0 0 0 0 50];
% X1 = simulate(mc,numSteps,'X0',x0);
% 
% figure;
% simplot(mc,X1);
end

% Find the most common sequences/transitions
T = find_patterns(VA);
%...  find specific interesting transition



% Single element transitions- get the lengths of each cell 
L = cellfun(@length,T);
%%extract cells of length 5 
iwant = T(L==2);


% Output variables:
out_markov.T =T; 
out_markov.VA = VA; % vector of flight identities
out_markov.FL_clust = out; % segregated flights
out_markov.day = day2use; % segregated flights





