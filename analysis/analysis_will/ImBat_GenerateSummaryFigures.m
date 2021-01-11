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
    
    col = hsv(7);
    for i = 1 : length(subFolders)
        load([subFolders(i).name, '/CellReg_files/ROI_Data/Saved_Data/Aligned_Data.mat'])
        
        for ii = 1:3
            % 1: FLight STATS
            combined.flights{i}.FL_consistancy{ii} = ImBat_FlightStats(FlightAlignedROI{ii});
        
            % 2: Calcium STATS
            combined.calcium{i}.ScoreMatrix{ii} = ImBat_BuildScoreMatrix(CombinedROI,FlightAlignedROI{ii});
        end
        
        
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
    
    figure();
    hold on;
    for ii = 1 : length(subFolders)
        % 1: Calcium plots
        for iii = 1:2%: length(combined.flights{i}.FL_consistancy);
            subplot(1,2,iii);
            hold on;
            ImBat_Tuning_Stability( combined.calcium{ii}.ScoreMatrix{iii},col(ii,:));
            grid on;
                xlim([0 6]);
                ylim([-0.2 1.05]);
        end
    end
    grid on
    
    
    
    %legend('Bat 1','Bat 1', 'Bat 1','Bat 1','Bat 1','Bat 1','Bat 1'
    
    % 2: Calcium plots
    
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

