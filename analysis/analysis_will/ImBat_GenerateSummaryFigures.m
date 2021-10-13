function [SumData] = ImBat_GenerateSummaryFigures;

% index into all processed folders ( desktop usually)
% For making all batch figures...

% WAL3
% 12/05/2020
extract_data = 1;
make_plots =1;
reWrite = 0; % overwrite data ( shoul be off...);

DIR = pwd;

files = dir(pwd);
files(ismember( {files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.

% Init data Variables:
FirstThreeClusters = [];
Clust2save = []; % for plotting clusters over time in session
AllFLights2save = []; % for plotting clusters over time in session
counter = 1;
if extract_data ==1;
    disp(' extracting and aligning data to flights');
    for i = 1 : length(subFolders)
        cd([subFolders(i).name, '/CellReg_files/ROI_Data']);
        
        % % CHECK, or possibly, Remove..
        if reWrite ==1;
            if exist('Saved_Data','dir')>6
                rmdir('Saved_Data','s');
                disp('Removed folder...');
            end
        end
        if exist('Saved_Data','dir')>6
            disp('Data already extracted, moving on');
        else
            tmp=dir(fullfile(pwd,'*.mat'));
            file_list={tmp(:).name};
            
            % load all mat files
            for ii = 1: length(file_list);
                load([file_list{ii}]);
            end
            
            % Concatonate/Cluster data across days
            [CombinedROI,ROI_Data] = ImBat_GroupCalcium(ROI_Data,cell_registered_struct,aligned_data_struct);
            
            % repair z bat flight data:
            ROI_Data = ImBat_RepairFlightData(ROI_Data);
            
            % now, we have a restricted set for ROI_Data, now cluster the flights:
            flightPaths = ImBat_GroupFlights(ROI_Data,'mtf',master_track_file,'dist',1.2);         % just the flights
            close all
            
            % align a few flights
            for ii = 1:3
                [FlightAlignedROI{ii}] = ImBat_Align_FC(CombinedROI,flightPaths,ii+1);
            end
            close all
            
            mkdir('Saved_Data')
            save('Saved_Data/Aligned_Data.mat','flightPaths','CombinedROI','FlightAlignedROI','-v7.3');
            clear flightPaths CombinedROI ROI_Data FlightAlignedROI
        end
        cd(DIR);
    end
    
end
% Now, you can load the data, and run batch proccessing scripts and make figures

if make_plots ==1
    disp('Perfoming Summary analysis');
    
    col = hsv(length(subFolders)+1);
    for i = 1 : length(subFolders)
        load([subFolders(i).name, '/CellReg_files/ROI_Data/Saved_Data/Aligned_Data.mat'])
        load([subFolders(i).name, '/CellReg_files/ROI_Data/ROI_Data.mat'])
        
        
        %% generate Null data:
        FlightAligned_FL_Sorted = FlightAlignedROI;  % random sort flights
        FlightAligned_ROI_Sorted = FlightAlignedROI; % random sort Ca ID
        FlightAlighed_Sorted = ImBat_SortByBehav(FlightAlignedROI); %  sort Ca by tuning
        
numItter = 5;

for iii = 1:3; % flight types
    % shuffle Day index
    h = FlightAlignedROI{iii}.CutCells_date;
    fakeInd = randperm((size(h,2)));
    FlightAligned_FL_Sorted{iii}.CutCells_date = FlightAlignedROI{iii}.CutCells_date(fakeInd);

for ivi = 1:numItter; % this will add an adendum of ROIs per day
    
for iv = 1:max(h) %  % shuffle ROI indexfor all ROIs
 idx2use = find(h==iv);
 h2a = size(FlightAlignedROI{iii}.C,1);
h2b = size(FlightAlignedROI{iii}.C,3);

fakeInd_1 = randperm(h2a);
%fakeInd_2 = randperm(h2b);

% just keep adding to the index
FlightAligned_ROI_Sorted{iii}.C(1+(size(FlightAlignedROI{iii}.C,1)*(ivi-1)):length(fakeInd_1)+(size(FlightAlignedROI{iii}.C,1)*(ivi-1)),:,idx2use) = FlightAlignedROI{iii}.C(fakeInd_1,:,idx2use);
FlightAligned_ROI_Sorted{iii}.C_raw(1+(size(FlightAlignedROI{iii}.C,1)*(ivi-1)):length(fakeInd_1)+(size(FlightAlignedROI{iii}.C,1)*(ivi-1)),:,idx2use) = FlightAlignedROI{iii}.C_raw(fakeInd_1,:,idx2use);
FlightAligned_ROI_Sorted{iii}.S(1+(size(FlightAlignedROI{iii}.C,1)*(ivi-1)):length(fakeInd_1)+(size(FlightAlignedROI{iii}.C,1)*(ivi-1)),:,idx2use) = FlightAlignedROI{iii}.S(fakeInd_1,:,idx2use);
end
end


clear h fakeInd h2 fakeInd_Ca  fakeInd_1 fakeInd_2
end
        
        
        %% Get aggregtaed stats
        
        for ii = 1:3
            % 1: FLight STATS
            combined.flights{i}.FL_consistancy{ii} = ImBat_FlightStats(FlightAlignedROI{ii});
            
            % 2: Calcium STATS
            combined.calcium{i}.ScoreMatrix{ii} = ImBat_BuildScoreMatrix(CombinedROI,FlightAlignedROI{ii});
            
            % 2b. Calcium STATS sorted by ROI tuning:
            combined.calcium_sorted{i}.ScoreMatrix{ii} = ImBat_BuildScoreMatrix(CombinedROI,FlightAlighed_Sorted{ii});
            % 2c. Calcium randominzing ROIs ( null):
            combined.calcium_Ca_shuffled{i}.ScoreMatrix{ii} = ImBat_BuildScoreMatrix(CombinedROI,FlightAligned_ROI_Sorted{ii});
            % 2d. Calcium randominzing Flights:
            combined.calcium_Fl_shuffled{i}.ScoreMatrix{ii} = ImBat_BuildScoreMatrix(CombinedROI,FlightAligned_FL_Sorted{ii});
            
            
            
        end
        
        close all
        % Get more(!) flight stats
        [Clust2save(:,:,i), AllFLights2save(:,:,i)] = ImBat_Session_FlightRate(flightPaths);
        
        
        FL_Stats = ImBat_Quantify_Flights(flightPaths);
        FirstThreeClusters = cat(1,FirstThreeClusters,FL_Stats.FirstThreeClusters);
        
        % combined.calcium = 0
        
    end
    
    %% Plotting... save in 'Processed/Figures'
    temp_FL_m = [];
    temp_FL_std = [];
    temp_FL_s = [];
    counter =1;
    figure();
    hold on;
    consol_flights.toPlot = [];
        consol_flights.pval_combined_data = [];

    for ii = 1 : length(subFolders)
        % 1: Flight plots
        for iii = 1:3%: length(combined.flights{i}.FL_consistancy);
            subplot(1,4,iii);
            ImBat_Plot_FlightStats(combined.flights{ii}.FL_consistancy{iii},col(ii,:));
            grid on;
            %consolidate flights
            if iii <4; % ==1;
            try
                temp_FL_m(counter,1:length(combined.flights{ii}.FL_consistancy{iii}.toPlot(1,:))) = combined.flights{ii}.FL_consistancy{iii}.toPlot(1,:);
                temp_FL_std(counter,1:length(combined.flights{ii}.FL_consistancy{iii}.toPlot(2,:))) = combined.flights{ii}.FL_consistancy{iii}.toPlot(2,:);
                temp_FL_s(counter,1:length(combined.flights{ii}.FL_consistancy{iii}.toPlot(3,:))) = combined.flights{ii}.FL_consistancy{iii}.toPlot(3,:);
                counter = counter+1;
            catch
                disp('error')
            end   
            end
        end
    end
    grid on
    
    temp_FL_m(temp_FL_m==0) = NaN;
    temp_FL_std(temp_FL_std==0) = NaN;
    ind2small = temp_FL_s<10; % remove really small indexes too.
    temp_FL_m(ind2small) = NaN;
    temp_FL_std(ind2small) = NaN;
    
    consol_flights.toPlot(1,:) = nanmedian(temp_FL_m);
    consol_flights.toPlot(2,:) = nanmedian(temp_FL_std);
    consol_flights.toPlot(3,:) = ones(1,size(mean(temp_FL_m),2));
    % plot all flights together:
    figure(); ImBat_Plot_FlightStats(consol_flights,col(ii,:));

    
    
    %% Plot calcium stuff:
    consol_calcium = [];
    consol_calcium_sorted = [];
    consol_calcium_Ca_shuffled = [];
    consol_calcium_Fl_shuffled = [];
    
    figure();
    hold on;
    for ii = 1 : length(subFolders)
        % 1: Calcium plots
        for iii = 1:3%: length(combined.flights{i}.FL_consistancy);
            hold on;
            ImBat_Tuning_Stability( combined.calcium{ii}.ScoreMatrix{iii},col(ii,:));
            grid on;
            xlim([0 11]);
            ylim([-0.4 1.05]);
            try
                consol_calcium = cat(2,consol_calcium,combined.calcium{ii}.ScoreMatrix{iii});
                consol_calcium_sorted = cat(2,consol_calcium_sorted,combined.calcium_sorted{ii}.ScoreMatrix{iii});
                
                consol_calcium_Ca_shuffled = cat(2,consol_calcium_Ca_shuffled,combined.calcium_Ca_shuffled{ii}.ScoreMatrix{iii});
                consol_calcium_Fl_shuffled = cat(2,consol_calcium_Fl_shuffled,combined.calcium_Fl_shuffled{ii}.ScoreMatrix{iii});
                
            catch
            end
        end
    end
    grid on
    
    figure();
    ImBat_Tuning_Stability(consol_calcium,col(ii,:))
    ImBat_Tuning_Stability(consol_calcium_sorted,col(ii-1,:))
    ImBat_Tuning_Stability(consol_calcium_Ca_shuffled,col(ii-2,:))
    ImBat_Tuning_Stability(consol_calcium_Fl_shuffled,col(ii-3,:))
    
    ylim([-0.4 1.05]);
    title('Tuning Stability to day 1');
    ylabel('Corr to day 1')
    
    % Display all stability:
    
    
    % Plot top FLights percentages/day:
    figure();
    FirstThreeClustersT = FirstThreeClusters;
    FirstThreeClustersT(FirstThreeClustersT==0) = NaN;
    plotSpread(FirstThreeClustersT,'showMM',4);
    title('Flightpath Sorted by session');
    
    %     % NOW, resort based on max for that day:
    FirstThreeClustersT2 = FirstThreeClusters(:,2:end);
    FirstThreeClustersT2 = sort(FirstThreeClustersT2,2,'descend');
    FirstThreeClustersT2 = cat(2,FirstThreeClustersT(:,1),FirstThreeClustersT2);
    %
    figure();
    FirstThreeClustersT2(FirstThreeClustersT2==0) = NaN;
    plotSpread(FirstThreeClustersT2,'showMM',4);
    title('Flightpath Sorted by day');
    
    %legend('Bat 1','Bat 1', 'Bat 1','Bat 1','Bat 1','Bat 1','Bat 1'
    
    %%  flights over session time:
    for i = 1: size(Clust2save,3);
        NormClust2save(:,:,i) =  squeeze(Clust2save(:,:,i))./squeeze(AllFLights2save(:,:,i));
    end
    NormClust2save2(1,:,:) = NormClust2save(1,:,:);
    NormClust2save2(2,:,:) = sum(NormClust2save(2:end,:,:),1);
    figure();
    hold on;
    
    for i = 1:size(NormClust2save2,1)
        clear adata
        adata = squeeze(NormClust2save2(i,1:50,:))';
        adata(isnan(adata)) = 0.1;
        L = size(adata,2);
        se = nanstd(adata)/2;%sqrt(length(adata));
        mn = nanmedian(adata);
        mn = smooth(mn,10)';
        se = smooth(se,10)';
        h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col(i,:)); alpha(0.5);
        plot(mn,'Color',col(i,:));
    end
    
    
end






%% Aggregate data/summary for exporting

% Basic summary data combined across animals
SumData.PlaceCells_per_day = []; % number of place cells/day
SumData.Cells_per_day = []; % number of cells/day
SumData.Flights_per_day = []; % number of total flights/day

%% Stability data aggregate
SumData.Stability.ROI_stability = consol_calcium; % number of place cells/day
% Sorted ROI
SumData.Stability.SortedROI_ROI_stability = consol_calcium_sorted;
% Shuffled ROIs
SumData.Stability.ShuffledROI_ROI_stability = consol_calcium_Ca_shuffled;
SumData.Stability.ShuffledROI_numItter = numItter; 
% Shuffled Flights
SumData.Stability.ShuffledFlights_ROI_stability = consol_calcium_Fl_shuffled;
SumData.Stability.ShuffledFlights_numItter = numItter; 


% number of flights...

% Now, get basic things about the flights:
[dat] = ImBat_Totals


% COMBINEING ANIMALS:
% 1. All place cell tracking ( how many can we track, and for how long)
% 2. Place cell stability ( first and last, blue/red dots) combining all animals
% 3. Place cell stability correlation
% 4. Top flight stability ( PCA of position)

