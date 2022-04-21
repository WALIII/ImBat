function [out_markov] = ImBat_New_Markov(flightPaths)

% Segregate Flights:
%day2use = 12;
Sim = 0;
% 
% [out] =  ImBat_SegTrajectories(ROI_Data{1,day2use}.Alignment.flightPaths.flights,ROI_Data{1,day2use}.Alignment.flightPaths.Location_time,'nclusters',nclust);
% close all
% 
% % Plot flights
% 
% figure(); 
% hold on;
% plot3(ROI_Data{1,day2use}.Alignment.flightPaths.flights(:,1),ROI_Data{1,day2use}.Alignment.flightPaths.flights(:,2),ROI_Data{1,day2use}.Alignment.flightPaths.flights(:,3),'k');% plot the flight trajectory in space
% 
% % plot all segregated trajectories
% hold on;
% for i = 1:size(flightPaths.flight_starts_idx,2)
% plot3(ROI_Data{1,day2use}.Alignment.flightPaths.flights(flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i),1),ROI_Data{1,day2use}.Alignment.flightPaths.flights(flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i),2),ROI_Data{1,day2use}.Alignment.flightPaths.flights(flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i),3),'r');% plot the flight trajectory in space
% end
% 

% Plot specific flights
figure();
hold on;
colorC = hsv(10);
try
for i = 2:10;
    for ii = 1: size(flightPaths.clusterIndex{i},1)
        idX = flightPaths.clusterIndex{i}(ii);
        plot3(flightPaths.trajectoriesContinous(1,flightPaths.flight_starts_idx(idX):flightPaths.flight_ends_idx(idX)),flightPaths.trajectoriesContinous(2,flightPaths.flight_starts_idx(idX):flightPaths.flight_ends_idx(idX)),flightPaths.trajectoriesContinous(3,flightPaths.flight_starts_idx(idX):flightPaths.flight_ends_idx(idX)),'Color',colorC(i,:));% plot the flight trajectory in space

    end
end
catch
end
grid on;







%ROI_Data{1,day2use}.Alignment.flightPaths.flights(flightPaths.flight_starts_idx(idX):flightPaths.flight_ends_idx(idX),2),ROI_Data{1,day2use}.Alignment.flightPaths.flights(flightPaths.flight_starts_idx(idX):flightPaths.flight_ends_idx(idX),3),'Color',colorC(i,:));% plot the flight trajectory in space


% Plot specific flights
figure();
hold on;
colorC = hsv(size(flightPaths.clusterIndex,2));
try
for i = 1:10%size(flightPaths.clusterIndex,2)
    subplot(5,2,i);
    hold on;
    for ii = 1: size(flightPaths.clusterIndex{i},1)
        idX = flightPaths.clusterIndex{i}(ii);
        plot3(flightPaths.trajectoriesContinous(1,flightPaths.flight_starts_idx(idX):flightPaths.flight_ends_idx(idX)),flightPaths.trajectoriesContinous(2,flightPaths.flight_starts_idx(idX):flightPaths.flight_ends_idx(idX)),flightPaths.trajectoriesContinous(3,flightPaths.flight_starts_idx(idX):flightPaths.flight_ends_idx(idX)),'Color',colorC(i,:));% plot the flight trajectory in space
    end
    title(['Flight ', num2str(i)]);
end
catch
end



% build transition matrix:
% sort based on times.. 


% Seq = zeros(1,size(flightPaths.flight_ends_idx,2));
% for i = 1:size(flightPaths.clusterIndex,2)
%     Seq(flightPaths.clusterIndex{i}) = i; %create single vector of flight IDs
%     Seq_time(flightPaths.clusterIndex{i}) = flightPaths.flight_starts_idx(flightPaths.clusterIndex{i});
% Col2use(i) = size(flightPaths.clusterIndex{i},1);
% end

Seq = flightPaths.id'; %create single vector of flight IDs
Seq_time = flightPaths.flight_starts_idx(:)';
for i = 1:size(flightPaths.clusterIndex,2)
Col2use(i) = size(flightPaths.clusterIndex{i},1);
end

[aa,ab] = sort(Seq_time);
Seq = Seq(ab);
Seq_time = Seq_time(ab);


% build sequence
%Seq2 = (1:size(flightPaths.flight_starts_idx,2));
VA(:,1) = Seq;
% VA(:,2) = Seq2;
 [u,~,n] = unique(VA,'rows');
 N = length(u);
 F = accumarray([n(1:end-1),n(2:end)],1,[N,N]);
 T = bsxfun(@rdivide,F,sum(F,2));
T(T<0.2) = 0;
 
 % Make simple markov model
 mc = dtmc(T);
 figure;
G = graphplot(mc,'ColorEdges',true);
G.MarkerSize = ceil(mat2gray(Col2use)*20)+1;


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

% % flight lengths 
% for i = 1:max(flightPaths.id)
% ID_len(i) = mean(flightPaths.flight_ends_idx(find(flightPaths.id==i))-flightPaths.flight_starts_idx(find(flightPaths.id==i)))/120;
% end

 % less robust:
% for i = 1:length(flightPaths.flight_ends_idx)
%     FL_len(i) = flightPaths.flight_ends_idx(i)-flightPaths.flight_starts_idx(i);
% end

% % more robust:
G = squeeze((flightPaths.vel(:,:,:)));
for i = 1:length(flightPaths.flight_ends_idx)
    glen = find(G(100:end,i)<0.3);
    FL_len(i) = glen(2)+100;
    clear glen
end

% Output variables:
out_markov.T =T; 
out_markov.VA = VA; % vector of flight identities
out_markov.out_sort = ab;

out_markov.FlightIDVector = Seq;
out_markov.FlightTimeVector = Seq_time;
% dont forget to sort!
% out_markov.ID_Len = ID_len; % length of each flight ID in time
out_markov.FL_len = FL_len(ab); %length of each individual flight in time




