function [ROI_Data] = ImBat_MultiDayAnalysis
% First-pass Multi-day analysis function 

% WAL3
% d190609

% User input: processing steps to run ( t/f)
ST1 = 1; % extract data
ST2 = 1; % multiday alignment
ST3 = 1; % Across day analysis
    ST3_1 = 1; % Max projection overlay 
    ST3_2 = 1; % Place Cell Overlay
    ST3_3 = 1; %
    ST3_5 = 1;

% local Directory:
LD =  'F:\Bat_search';
%'/Users/ARGO/Documents/DATA/Bat_Data';

% Cell array containing strings of the days you want to look at
days = {'190528','190529','190530'};

% Run in the 'processed' folder containing the days
DIR = pwd;

% check into DIRs
cd(LD);
cd(DIR);

%% Go through all days and grab the data
if ST1 == 1;
for i = 1: size(days,2)
    % Save the date:
    ROI_Data{i}.date = days{i};
    
    % create index into correct subfolder
    folder = [days{i},'/','Zack_',days{i},'_fly-1_extraction/processed'];
    folder_rest = [days{i},'/','Zack_',days{i},'_rest1-1_extraction/processed'];
    
    disp(['entering into folder for day:  ', days{i}]);
    cd(folder) % index into folder 
    Alignment = load('Alignment.mat'); % get alignment data
    ROIs = load('Motion_corrected_Data_DS_results.mat'); % Get ROI data
   
    % Consolidate ROI and Flight data
    ROI_Data{i}.Alignment = Alignment;
    ROI_Data{i}.ROIs = ROIs;
    clear Alignment ROIs

    % Make Max projections ( Flights )
    disp( 'Getting Max Projection for Flights...');
    load('Motion_corrected_Data_DS','Y','Ysiz');
    [MaxProj_flights, ~] = ImBat_Dff(Y);
    ROI_Data{i}.MaxProj_flights = MaxProj_flights;
    ROI_Data{i}.Ysiz_full = Ysiz;
    clear Y Ysiz MaxProj_flights;
    
    % Make Max projections ( Rest )
    cd(DIR); cd(folder_rest);
    load('Motion_corrected_Data_DS', 'Y','Ysiz');
    disp( 'Getting Max Projection for rest...');
    [MaxProj_rest, ~] = ImBat_Dff(Y);
    ROI_Data{i}.MaxProj_rest = MaxProj_rest;
    ROI_Data{i}.Ysiz_DS = Ysiz;
    clear Y Ysiz MaxProj_rest;

        cd(DIR)
end

% Save data locally somewhere...
cd(LD);
for i = 1:size(ROI_Data,2);
    if i ==1;
        filename = ['ROI_Data_',days{1}]
    else
        filename = [filename,'_',days{i}];
    end
end
    
save([filename,'.mat'],'ROI_Data')

end

%% Cell Sort to get 'similar ROIs'
if ST2 ==1;
    
end


%% Plot place cell-ness of the top most consistant cells
if ST3 == 1
 % Basic alignment within day:
 figure();
 for i = 1:3
     subplot(1,3,i);
 A = ROI_Data{i}.MaxProj_rest;
 B = ROI_Data{i}.MaxProj_flights;
[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = 0.009;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;
A = imregister(A, B, 'rigid', optimizer, metric);

temp = cat(3,A,B);
D_reg{i} = max(temp,[],3);
clear temp
[RGB1 RGB2] = CaBMI_XMASS(A,B,A);
imagesc(squeeze(RGB1))
 end
 
 clear A B C
 % Across day alignment ( rest)
 figure(); 
 
 A = ROI_Data{1}.MaxProj_rest;
 B = ROI_Data{2}.MaxProj_rest;
 C = ROI_Data{3}.MaxProj_rest;
 
[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = 0.009;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;
A = imregister(A, B, 'rigid', optimizer, metric);
C = imregister(C, B, 'rigid', optimizer, metric);
[RGB1 RGB2] = CaBMI_XMASS(A,B,C);
imagesc(squeeze(RGB1))


 clear A B C
 % Across day alignment ( flights)
 figure(); 
 
 A = ROI_Data{1}.MaxProj_flights;
 B = ROI_Data{2}.MaxProj_flights;
 C = ROI_Data{3}.MaxProj_flights;
 
A = imregister(A, B, 'rigid', optimizer, metric);
C = imregister(C, B, 'rigid', optimizer, metric);
[RGB1 RGB2] = CaBMI_XMASS(A,B,C);
imagesc(squeeze(RGB1))

clear AB
% across day alignemnt( all)
 figure(); 
 
 A = D_reg{1};
 B = D_reg{2};
 C = D_reg{3};
 
A = imregister(A, B, 'rigid', optimizer, metric);
C = imregister(C, B, 'rigid', optimizer, metric);
[RGB1 RGB2] = CaBMI_XMASS(A,B,C);
imagesc(squeeze(RGB1))



% Condition Registration by saving data in the proper format
mkdir('CellReg_files');
cd('CellReg_files')

for i = 1:size(days,2)
    
    % get size of ROI field
    A2 = full(ROI_Data{i}.ROIs.results.A);
    for ii = 1:size(A2,2);
        A3(ii,:,:) = reshape(A2(:,ii),ROI_Data{i}.Ysiz_DS(1),ROI_Data{i}.Ysiz_DS(2)); 
    end
         A2 = A3(1:100,:,:); % ROI masks in CellReg format...

    % save data:
    save([ROI_Data{i}.date,'.mat'],'A2');
    clear A2 A3
    
end


mkdir('ROI_Data');
save('ROI_Data/ROI_Data','ROI_Data');
% Register Cells by ROI masks:
CellReg

% Save indexes for similar ROIs, and load into ROI_Data ( for across-day
% indexing)



end







