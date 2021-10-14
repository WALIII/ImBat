
function [R,Map2save, occupancyMap2save] = ImBat_Compare2D_ratemeaps_LightDark(FlightAlignedROI,out_data, roi2plot,cluster);


Flights = [];
Spikes =  [];
bound2use = 1:1400;

if size(FlightAlignedROI,2)>2
    disp('WARNING: make sure the ~Unique~ cluster has been added to FlightAlignedROI');
end
% cluster = 1;
condition2use = 1:3; % light or dark

for iiii = 1:length(condition2use);
    
    for ii = roi2plot;%1:size(CutCells,1);
        
        
        cell2use = ii;
        %ind2use = find(FlightAlignedROI{cluster}.CutCells_date == days2use(iiii));
        if iiii == 1; % first light condition
            ind2use = out_data.Light_cluster{cluster}.FirstLight;
        elseif iiii == 2; % first dark condition
            ind2use = out_data.Light_cluster{cluster}.Darkness;
        elseif iiii == 3; % second light condition
            ind2use = out_data.Light_cluster{cluster}.LastLight;
        end
        
        for ix = 1:size(ind2use,1); % correct the index
            ind_true(ix) = find(FlightAlignedROI{cluster}.cluster_idX== ind2use(ix));
        end
        ind2use = ind_true;
        % get flight data,
        ClustFlight = FlightAlignedROI{cluster}.ClustFlight_withPads(:,:,ind2use);
        
        % upsample the calcium
        CutCells = FlightAlignedROI{cluster}.S(:,:,ind2use);
        
        for i = 1: size(CutCells,3);
            
            
            colormap(hot);
            trial2use = i;
            exampFlight = ClustFlight(bound2use,:,trial2use);
            exampCell = CutCells(cell2use,:,trial2use);
            exampCell = interp(exampCell,4);
            exampCell = exampCell(1:size(exampFlight,1));
            
            if sum(exampCell)>0.00001
                Flights = cat(1,Flights, exampFlight);
                Spikes = cat(1,Spikes, exampCell');
                
            else
                % disp('not enought spikes');
            end
            
        end
    end
    
    
    % binarie spikes
    %     Spikes = zscore(Spikes);
    %     Spikes(Spikes>1) = 1;
    %     Spikes(Spikes<1) = 0;
    %     Spk = find(Spikes ==1);
    %     val = ones(length(Spk),1);
    %%%  Count Spikes
    Spikes = zscore(Spikes);
    Spikes = round(Spikes);
    [Spk] = find(Spikes>1);
    val = Spikes(Spk);
    
    
    try
        [NormRateMat, occupancyMap] =  ImBat_2dHeatMap(Flights',Spk',val);
        figure(); imagesc(((NormRateMat')))
        Map2save(:,:,iiii) = NormRateMat';
        occupancyMap2save(:,:,iiii) = occupancyMap';
    catch
    end
    clear Spikes Flights exampCell CutCells exampFlight ind2use
    % re-initialize vars..
    Flights = [];
    Spikes = [];
    
end

% compare across days

% calculate the first, within day:

for iiii = 1:2
    
    for ii = roi2plot;%1:size(CutCells,1);
        
        
        cell2use = ii;
        if iiii == 1; % first light condition
            ind2use = out_data.Light_cluster{cluster}.FirstLight;
        elseif iiii == 2; % first dark condition
            ind2use = out_data.Light_cluster{cluster}.Darkness;
        elseif iiii == 3; % second light condition
            ind2use = out_data.Light_cluster{cluster}.LastLight;
        end
        for ix = 1:size(ind2use,1); % correct the index
            ind_true(ix) = find(FlightAlignedROI{cluster}.cluster_idX== ind2use(ix));
        end
        ind2use = ind_true;
        
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
            
            colormap(hot);
            trial2use = i;
            exampFlight = ClustFlight(bound2use,:,trial2use);
            exampCell = CutCells(cell2use,:,trial2use);
            exampCell = interp(exampCell,4);
            exampCell = exampCell(1:size(exampFlight,1));
            
            if sum(exampCell)>0.00001
                Flights = cat(1,Flights, exampFlight);
                Spikes = cat(1,Spikes, exampCell');
            else
                %      disp('not enought spikes');
            end
            
        end
    end
    
    
    %%%  binarize spikes
    %     Spikes = zscore(Spikes);
    %     Spikes(Spikes>1) = 1;
    %     Spikes(Spikes<1) = 0;
    %     Spk = find(Spikes ==1);
    %     val = ones(length(Spk),1);
    % count Spikes
    Spikes = zscore(Spikes);
    Spikes = round(Spikes);
    [Spk] = find(Spikes>1);
    val = Spikes(Spk);
    
    try
        [NormRateMat, occupancyMap] =  ImBat_2dHeatMap(Flights',Spk',val);
        
        figure(); imagesc(((NormRateMat')))
        Map2save_EO(:,:,iiii) = NormRateMat';
        clear Spikes Flights exampCell CutCells exampFlight ind2use
        
    catch
    end
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

