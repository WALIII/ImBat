function [S2save,F2save,other] =  ImBat_Spikes(FlightAlignedROI);
% Align spikes and calcium to the same timescale by upsampling Ca data

% Used for various functions
NANPAD = 0;

% Default is first cluster


for ii = 1:size(FlightAlignedROI{1}.S,1) % take the first 100 ROIs
    Flights = [];
    Spikes =  [];
    T_forward = [];
    T_reverse = [];
    T_reward = [];
    T_spd = [];
    T_dist = [];
    
    hold on;
    bound2use = 1:1400;
    
    cell2use = ii;
    
    % get flight data,
    ClustFlight = FlightAlignedROI{1}.ClustFlight_withPads;
    
    % upsample the calcium
    CutCells = FlightAlignedROI{1}.S;
    
    for i = 1: size(CutCells,3);
        
        bound2use = 1:1400;
        trial2use = i;
        exampFlight = ClustFlight(bound2use,:,trial2use);
        exampCell = CutCells(cell2use,:,trial2use);
        d1 = exampCell;
        exampCell = interp(exampCell,4);
        exampCell = exampCell(1:size(exampFlight,1));
        
        % get rid of edges
        exampCell(1:500) = [];
        exampFlight(1:500,:) = [];
        
        exampCell(end-300:end) = [];
        exampFlight(end-300:end,:) = [];
        time2export_forward = 1:length(exampFlight);
        time2export_reverse = fliplr(1:length(exampFlight));
        reward = zeros(1,length(exampFlight));
        reward(525:end-25) = 1;
        
        if NANPAD ==1;
            if sum(d1)==0
                exampCell = zeros(1,length(exampCell));
                exampCell(exampCell ==0) = NaN;
            end
        end
        
        % Calculate speed:
x1 = exampFlight(2:end,1);
x2 = exampFlight(1:end-1,1);
y1 = exampFlight(2:end,2);
y2 = exampFlight(1:end-1,2);
z1 = exampFlight(2:end,3);
z2 = exampFlight(1:end-1,3);

S=sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);
S(S>40) = 0;
spd = smooth(squeeze(S),20)';
dist = cumsum(spd);
        
        
% TO DO Direct distance to the reward:
        
        
        Flights = cat(1,Flights, exampFlight);
        Spikes = cat(1,Spikes, exampCell');
        T_forward = cat(1,T_forward,time2export_forward');
        T_reverse =  cat(1,T_reverse,time2export_reverse');
        T_reward = cat(1,T_reward,reward');
        T_spd = cat(1,T_spd,spd',spd(end));
        T_dist =  cat(1,T_dist,dist',dist(end));
     
    end
    % binarie spikes
    Spikes =  smooth(zscore(Spikes),10);
    % Spikes(Spikes>1) = 1;
     Spikes(Spikes<0) = 0;
    
    
    S2save(:,ii) = Spikes;
    F2save = Flights;
    other.T_forward = T_forward;
    other.T_reverse = T_forward;
    other.T_reward = T_reward;
    other.T_dist = T_dist;
    other.T_spd = T_spd;
end
