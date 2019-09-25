function ImBat_Markov(ROI_Data)

% Segregate Flights:

[out] =  ImBat_SegTrajectories(ROI_Data{1,12}.Alignment.out.flights,ROI_Data{1,12}.Alignment.out.Location_time);


% Plot flights

figure(); 
hold on;
plot3(ROI_Data{1,12}.Alignment.out.flights(:,1),ROI_Data{1,12}.Alignment.out.flights(:,2),ROI_Data{1,12}.Alignment.out.flights(:,3),'k');% plot the flight trajectory in space

% plot all segregated trajectories
hold on;
for i = 1:67
plot3(ROI_Data{1,12}.Alignment.out.flights(out.flight_starts_indx(i):out.flight_ends_indx(i),1),ROI_Data{1,12}.Alignment.out.flights(out.flight_starts_indx(i):out.flight_ends_indx(i),2),ROI_Data{1,12}.Alignment.out.flights(out.flight_starts_indx(i):out.flight_ends_indx(i),3),'r');% plot the flight trajectory in space
end


% Plot specific flights
figure();
hold on;
colorC = hsv(size(out.ClusterIndex,2));
for i = 1:size(out.ClusterIndex,2)
    for ii = 1: size(out.ClusterIndex{i},2)
        idX = out.ClusterIndex{i}(ii);
        plot3(ROI_Data{1,12}.Alignment.out.flights(out.flight_starts_indx(idX):out.flight_ends_indx(idX),1),ROI_Data{1,12}.Alignment.out.flights(out.flight_starts_indx(idX):out.flight_ends_indx(idX),2),ROI_Data{1,12}.Alignment.out.flights(out.flight_starts_indx(idX):out.flight_ends_indx(idX),3),'Color',colorC(i,:));% plot the flight trajectory in space
    end
end




% plot all segregated trajectories
hold on;
for i = 1:67
plot3(ROI_Data{1,12}.Alignment.out.flights(out.flight_starts_indx(i):out.flight_ends_indx(i),1),ROI_Data{1,12}.Alignment.out.flights(out.flight_starts_indx(i):out.flight_ends_indx(i),2),ROI_Data{1,12}.Alignment.out.flights(out.flight_starts_indx(i):out.flight_ends_indx(i),3),'r');% plot the flight trajectory in space
end




% Plot specific flights
figure();
hold on;
colorC = hsv(size(out.ClusterIndex,2));
for i = 1:size(out.ClusterIndex,2)
    subplot(5,2,i);
    hold on;
    for ii = 1: size(out.ClusterIndex{i},2)
        idX = out.ClusterIndex{i}(ii);
        plot3(ROI_Data{1,12}.Alignment.out.flights(out.flight_starts_indx(idX):out.flight_ends_indx(idX),1),ROI_Data{1,12}.Alignment.out.flights(out.flight_starts_indx(idX):out.flight_ends_indx(idX),2),ROI_Data{1,12}.Alignment.out.flights(out.flight_starts_indx(idX):out.flight_ends_indx(idX),3),'Color',colorC(i,:));% plot the flight trajectory in space
    end
    title(['Flight ', num2str(i)]);
end




% build transition matrix:
Seq = zeros(1,size(out.flight_ends_indx,2));
for i = 1:size(out.ClusterIndex,2)
    Seq(out.ClusterIndex{i}) = i;
end


% build sequence
Seq2 = (1:67);
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
numSteps = 50;
X = simulate(mc,numSteps);
figure(); simplot(mc,X,'FrameRate',0.5,'Type','graph');



x0 = [0 0 0 0 0 50 0 0];
X1 = simulate(mc,numSteps,'X0',x0);

figure;
simplot(mc,X1);