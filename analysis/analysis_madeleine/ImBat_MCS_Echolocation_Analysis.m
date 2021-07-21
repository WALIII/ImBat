function [gen_echo_stats] = ImBat_MCS_Echolocation_Analysis(flightPaths)

% 1. Are there more or fewer echolocations on the stereotyped flights?
% 2. Does echolocation frequency change as a f(x) of session duration?
% 3. Does echolocation occurrence change as a f(x) of training duration?
% 4. Are echolocations on stereotyped flight in the same locations on every trial?


% 1. For every echolocation, what is the probability it has ocurred during
% a stereotyped flight?

% Create randomly spread distributio of echolocation calls (same #)
% Compare the probability of an echoloication occurring on a stereotyped flight to the probability of a FAKE echolocation occurring on a stereotyped flight

vectorized_JT = squeeze(flightPaths34.JT);
VJT = reshape(vectorized_JT,[],1);
VJT(isnan(VJT))=0;
no_nan_VJT = VJT(VJT~=0);
figure(); hold on; title("Real Data Echolocation Peaks All Flights"); plot(VJT);

ST = 0; UST = 0;
for j=1:size(VJT,1)
    if VJT(j) == 0.01
        flight_idx = ceil(j/1182);
        if flightPaths34.id(flight_idx) > 1
            ST = ST+1;
        else
            UST = UST+1;
        end
    end
end

PST = ST/size(no_nan_VJT,1);
PUST = UST/size(no_nan_VJT,1);

% Create randomly uniform echolocation vector of same size as real data
for i=1:1000
    tt = find(VJT==0.01);
    rs = 1;
    re = max(tt);
    r = round((re-rs).*rand(size(no_nan_VJT,1),1) + rs);
    rand_VJT = zeros(size(VJT,1),1);
    rand_VJT(r)=0.01;
    %figure(); plot(rand_VJT); hold on; title("SIM Data Echolocation Peaks All Flights");

    ST = 0; UST = 0;
    for j=1:size(rand_VJT,1)
        if rand_VJT(j) == 0.01
            flight_idx = ceil(j/1182);
            if flightPaths34.id(flight_idx) > 1
                ST = ST+1;
            else
                UST = UST+1;
            end
        end
    end

    SIM_PST(i) = ST/size(no_nan_VJT,1);
    SIM_PUST(i) = UST/size(no_nan_VJT,1);
end

MSIM_PST = mean(SIM_PST);
MSIM_PUST = mean(SIM_PUST);

% 2. Does echolocation occurance change as a f(x) of session duration?
flight_idx = [];
for j=1:size(VJT,1)
    if VJT(j) == 0.01
        flight_idx = [flight_idx,ceil(j/1182)];
    end
end

PST = ST/size(no_nan_VJT,1);
PUST = UST/size(no_nan_VJT,1);

% Where in the flight do they tend to echolocate?
for i=1:size(flightPaths34.id,1)
    figure(); hold on;
    plot3(flightPaths34.pos(1,:,i),flightPaths34.pos(2,:,i),flightPaths34.pos(3,:,i),'Color','b');
    subset = squeeze(flightPaths34.echos(:,:,i));
    try
        scatter3(flightPaths34.pos(1,find(subset==0.01),i),flightPaths34.pos(2,find(subset==0.01),i),flightPaths34.pos(3,find(subset==0.01),i),[],'*r');
    catch
        disp("This flightpath has no echolocations");
    end
end

%% Calculate curvature and correlate curvature to echolocations
for i=40:size(flightPaths34.id,1)-2
    X = flightPaths34.pos(:,:,i)';
    Xidx = find(~isnan(sum(X')));
    X_trunc = X(Xidx(1):Xidx(end-1),:);
    %figure(); hold on; 
    %scatter3(X_trunc(:,1),X_trunc(:,2),X_trunc(:,3),[],'ob');
    X_trunc=X_trunc(50:end-50,:);
    % Interpolate across NaN values
    nanx = isnan(X_trunc);
    nan_t = 1:numel(X_trunc);
    X_trunc(nanx) = interp1(nan_t(~nanx), X_trunc(~nanx), nan_t(nanx));
    to_delete = [];
    for j=1:size(X_trunc,1)-1
        if X_trunc(j,:) == X_trunc(j+1,:)
            to_delete = [to_delete,j];
        end
    end
    X_trunc(to_delete,:) = [];
    % Resample uniformly?
    interparc_t = size(flightPaths34.echos,2);
    X_trunc_i = interparc(interparc_t,squeeze(X_trunc(:,1)),squeeze(X_trunc(:,2)),squeeze(X_trunc(:,3)));
    %scatter3(X_trunc_i(:,1),X_trunc_i(:,2),X_trunc_i(:,3),[],'*r');
    
    % Get curvature measurements
    [L,R,K] = curvature(X_trunc_i);
    %figure();
    %plot(L,R); xlabel L; ylabel R
    %title('Curvature radius vs. cumulative curve length')
    figure();
    h = plot3(X_trunc_i(:,1),X_trunc_i(:,2),X_trunc_i(:,3)); 
    grid on; axis equal;
    set(h,'marker','.'); xlabel x; ylabel y; zlabel z
    title('3D curve with curvature vectors')
    hold on
    quiver3(X_trunc_i(:,1),X_trunc_i(:,2),X_trunc_i(:,3),K(:,1),K(:,2),K(:,3),10);
    hold off
    %disp(round(sum(abs(K(~isnan(K))))));
    % Plot CDF of curvature
    running_sum_R=[];
    for j=1:size(R,1)
        if j==1
            if isnan(R(j))
                running_sum_R(j) = 0;
            else
                running_sum_R(j) = R(1);
            end
        elseif ~isnan(R(j))
            running_sum_R(j) = running_sum_R(j-1)+R(j);
        elseif isnan(R(j))
            running_sum_R(j) = running_sum_R(j-1);
        end
    end
    figure(); hold on; %plot(running_sum_R);
    plot(R);
    echo_vect = squeeze(flightPaths34.echos(1,:,i));
    echo_vect(echo_vect==0)=NaN;
    echo_vect(echo_vect==0.01)=290000;
    stem(echo_vect,'r');
    title('Curvature radius (greater means more curved) with echolocations');      
end
    


end