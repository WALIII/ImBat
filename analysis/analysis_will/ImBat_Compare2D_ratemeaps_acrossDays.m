
function [R,Map2save,occupancyMap2save] = ImBat_Compare2D_ratemeaps_acrossDays(FlightAlignedROI,roi2plot,cluster);


biner = 1; % binerized data? 1 == yes, 2 == no
Flights = [];
Spikes =  [];

% cluster = 1;
days2use = unique(FlightAlignedROI{cluster}.CutCells_date);

NormRateMat = NaN(50,50);
occupancyMap = zeros(50,50);

for iiii = 1:length(days2use);
    
    for ii = roi2plot;%1:size(CutCells,1);
        
        bound2use = 1:1400;
        
        cell2use = ii;
        ind2use = find(FlightAlignedROI{cluster}.CutCells_date == days2use(iiii));
        
        % get flight data,
        ClustFlight = FlightAlignedROI{cluster}.ClustFlight_withPads(:,:,ind2use);
        
        % upsample the calcium
        CutCells = FlightAlignedROI{cluster}.S(:,:,ind2use);
        
        for i = 1: size(CutCells,3);
            
            bound2use = 1:1400;
            colormap(hot);
            trial2use = i;
            exampFlight = ClustFlight(bound2use,:,trial2use);
            exampCell = CutCells(cell2use,:,trial2use);
            exampCell = interp(exampCell,4);
            exampCell = exampCell(1:size(exampFlight,1));
            % get rid of edges
            exampCell(1:500) = 0;
            exampCell(end-300:end) = 0;
            Flights = cat(1,Flights, exampFlight);
            Spikes = cat(1,Spikes, exampCell');
            
        end
    end
    
    
    if biner ==1; % binerize spikes
        Spikes = zscore(Spikes);
        Spikes(Spikes>1) = 1;
        Spikes(Spikes<1) = 0;
        Spk = find(Spikes ==1);
        val = ones(length(Spk));
        
    else
        Spikes = zscore(Spikes);
        Spikes = round(Spikes);
        [Spk] = find(Spikes>1);
        val = Spikes(Spk);
    end
    try
        [NormRateMat, occupancyMap] =  ImBat_2dHeatMap(Flights',Spk',val);
    catch
        disp('');
    end
    %figure(); imagesc(((NormRateMat')))
    Map2save(:,:,iiii) = NormRateMat';
    occupancyMap2save(:,:,iiii) = occupancyMap';
    NormRateMat = NaN(50,50);
    
    clear Spikes Flights exampCell CutCells exampFlight ind2use
    % re-initialize vars..
    Flights = [];
    Spikes = [];
    
end

% compare across days

% calculate the first, within day:

for iiii = 1:2
    
    for ii = roi2plot;%1:size(CutCells,1);
        
        bound2use = 1:1400;
        
        cell2use = ii;
        ind2use = find(FlightAlignedROI{cluster}.CutCells_date == days2use(1));
        
        if iiii == 1;
            ind2use = ind2use(1:2:length(ind2use));
        elseif iiii ==2;
            ind2use = ind2use(2:2:length(ind2use));
        end
        
        
        % get flight data,
        ClustFlight = FlightAlignedROI{cluster}.ClustFlight_withPads(:,:,ind2use);
        
        % upsample the calcium
        CutCells = FlightAlignedROI{cluster}.S(:,:,ind2use);
        
        for i = 1: size(CutCells,3);
            
            bound2use = 1:1400;
            colormap(hot);
            trial2use = i;
            exampFlight = ClustFlight(bound2use,:,trial2use);
            exampCell = CutCells(cell2use,:,trial2use);
            exampCell = interp(exampCell,4);
            exampCell = exampCell(1:size(exampFlight,1));
            % get rid of edges
            exampCell(1:500) = 0;
            exampCell(end-300:end) = 0;
            Flights = cat(1,Flights, exampFlight);
            Spikes = cat(1,Spikes, exampCell');
            
        end
    end
    
    
    if biner ==1; % binerize spikes
        Spikes = zscore(Spikes);
        Spikes(Spikes>1) = 1;
        Spikes(Spikes<1) = 0;
        Spk = find(Spikes ==1);
        val = ones(length(Spk));
    else
        Spikes = zscore(Spikes);
        Spikes = round(Spikes);
        [Spk] = find(Spikes>1);
        val = Spikes(Spk);
    end
    
    NormRateMat =  ImBat_2dHeatMap(Flights',Spk',val);
   % figure(); imagesc(((NormRateMat')))
    
    Map2save_EO(:,:,iiii) = NormRateMat';
    NormRateMat = NaN(50,50);
    clear Spikes Flights exampCell CutCells exampFlight ind2use
    % re-initialize vars..
    Flights = [];
    Spikes = [];
    
end



for i = 1:size(Map2save,3)
    R(i) = corr2(squeeze(Map2save(:,:,1)),squeeze(Map2save(:,:,i)));
end

for i = 1
    R(1) = corr2(squeeze(Map2save_EO(:,:,1)),squeeze(Map2save_EO(:,:,2)));
end

% figure(); plot(R);

