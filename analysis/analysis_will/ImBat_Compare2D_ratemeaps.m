
function ImBat_Compare2D_ratemeaps(FlightAlignedROI,roi2plot);


Flights = [];
Spikes =  [];
iii = 1;

for ii = roi2plot;%1:size(CutCells,1);
hold on;
bound2use = 1:1400;

cell2use = ii;


% get flight data, 
ClustFlight = FlightAlignedROI{iii}.ClustFlight_withPads;

% upsample the calcium
CutCells = FlightAlignedROI{iii}.S;
%CutCells = FlightAlignedROI{iii}.S;


for i = 1: size(CutCells,3);
    

    bound2use = 1:1400;
    colormap(hot);
trial2use = i;
exampFlight = ClustFlight(bound2use,:,trial2use);
exampCell = CutCells(cell2use,:,trial2use);
exampCell = interp(exampCell,4);
exampCell = exampCell(1:size(exampFlight,1));

Flights = cat(1,Flights, exampFlight);
Spikes = cat(1,Spikes, exampCell');

end
end

% binarie spikes
% Spikes = zscore(Spikes);
% Spikes(Spikes>1) = 1;
% Spikes(Spikes<1) = 0;
  %%%  Count Spikes
    Spikes = zscore(Spikes);
    Spikes = round(Spikes);
    [Spk] = find(Spikes>1);
    val = Spikes(Spk);
    

Spk = find(Spikes ==1);
NormRateMat =  ImBat_2dHeatMap(Flights',Spk',val);
figure(); imagesc(((NormRateMat')))
