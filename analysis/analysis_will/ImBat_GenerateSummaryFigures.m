function ImBat_GenerateSummaryFigures;

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
        
        for ii = 1:3
            % 1: FLight STATS
            combined.flights{i}.FL_consistancy{ii} = ImBat_FlightStats(FlightAlignedROI{ii});
            
            % 2: Calcium STATS
            combined.calcium{i}.ScoreMatrix{ii} = ImBat_BuildScoreMatrix(CombinedROI,FlightAlignedROI{ii});
            
            % 3. Calcium STATS sorted by behav:
            FlightAlighed_Sorted = ImBat_SortByBehav(FlightAlignedROI);
           combined.calcium_sorted{i}.ScoreMatrix{ii} = ImBat_BuildScoreMatrix(CombinedROI,FlightAlighed_Sorted{ii});

        end
        
        close all
        % Get more(!) flight stats
        [Clust2save(:,:,i), AllFLights2save(:,:,i)] = ImBat_Session_FlightRate(flightPaths);
        
        
        FL_Stats = ImBat_Quantify_Flights(flightPaths,ROI_Data);
        FirstThreeClusters = cat(1,FirstThreeClusters,FL_Stats.FirstThreeClusters);
        
        % combined.calcium = 0
        
    end
    
    %% Plotting... save in 'Processed/Figures'
    figure();
    hold on;
    for ii = 1 : length(subFolders)
        % 1: Flight plots
        for iii = 1:3%: length(combined.flights{i}.FL_consistancy);
            subplot(1,4,iii);
            ImBat_Plot_FlightStats(combined.flights{ii}.FL_consistancy{iii},col(ii,:));
            grid on;
            %xlim([0 6]);
            
        end
    end
    grid on
    
    
    %% Plot calcium stuff:
    consol_calcium = [];
    consol_calcium_sorted = [];
    
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

            catch
            end
        end
    end
    grid on
    
    figure();
    ImBat_Tuning_Stability(consol_calcium,col(ii,:))
    ImBat_Tuning_Stability(consol_calcium_sorted,col(ii-1,:))

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


% 1. Plot Place cells as PNG
% 2. Plot flights
% 3. Plot heat map on flights for top 3 flights
% 4. Markov model (figure)
% 5. Higher order markov model (figure)
% 3. Get summary stats:
%     a.
%     b.
%     c,

% clear data, continue loop



% Aggregate data/summary

% PER ANIMAL: ( plot seperately, same figure though)
% number of place cells/day
% number of cells/day
% number of total flights/day
% number of flights...


% COMBINEING ANIMALS:
% 1. All place cell tracking ( how many can we track, and for how long)
% 2. Place cell stability ( first and last, blue/red dots) combining all animals
% 3. Place cell stability correlation
% 4. Top flight stability ( PCA of position)

