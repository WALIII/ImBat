extractedPlaceCells = dir('*ExtractedPlaceCells*');
load(extractedPlaceCells.name);
cd([pwd filesep 'Spatial_information']);

nClust =size(placeCells.pp_cells,1)-1; %number of clusters to look at excluding the first cluster
%initialize cells
cells_prePlacePost = cell(nClust,1);
cells_prePlace = cell(nClust,1);
cells_prePost = cell(nClust,1);
cells_placePost = cell(nClust,1);
cells_pre = cell(nClust,1);
cells_place = cell(nClust,1);
cells_post = cell(nClust,1);

%check if directories already exist and make if not
if ~exist([pwd filesep 'cells_prePlacePost'])
    mkdir([pwd filesep 'cells_prePlacePost'])
end
if ~exist([pwd filesep 'cells_prePlacePost' filesep 'prePlacePost'])
    mkdir([pwd filesep 'cells_prePlacePost' filesep 'prePlacePost'])
end
if ~exist([pwd filesep 'cells_prePlacePost' filesep 'prePlace'])
    mkdir([pwd filesep 'cells_prePlacePost' filesep 'prePlace'])
end
if ~exist([pwd filesep 'cells_prePlacePost' filesep 'prePost'])
    mkdir([pwd filesep 'cells_prePlacePost' filesep 'prePost'])
end
if ~exist([pwd filesep 'cells_prePlacePost' filesep 'placePost'])
    mkdir([pwd filesep 'cells_prePlacePost' filesep 'placePost'])
end
if ~exist([pwd filesep 'cells_prePlacePost' filesep 'pre'])
    mkdir([pwd filesep 'cells_prePlacePost' filesep 'pre'])
end
if ~exist([pwd filesep 'cells_prePlacePost' filesep 'place'])
    mkdir([pwd filesep 'cells_prePlacePost' filesep 'place'])
end
if ~exist([pwd filesep 'cells_prePlacePost' filesep 'post'])
    mkdir([pwd filesep 'cells_prePlacePost' filesep 'post'])
end

for clust_i = 1:nClust
    %pull out the putative pre,place,post cells from logical
    pp_cells = find(placeCells.pp_cells(clust_i+1,:));
    ppre_cells = find(placeCells.ppre_cells(clust_i+1,:));
    ppost_cells = find(placeCells.ppost_cells(clust_i+1,:));
    
    %find the different combos of pre, place, and post cells
    cells_prePlacePost{clust_i,:} = intersect(intersect(pp_cells,ppre_cells),ppost_cells);
    for cell_i = 1:length(cells_prePlacePost{clust_i,:})
        sourceFile = dir(['*_ROI_' num2str(cells_prePlacePost{clust_i}(cell_i)) '_cluster_' num2str(clust_i+1) '.tif']);
        if ~isempty(sourceFile)
            movefile(sourceFile.name,[pwd filesep 'cells_prePlacePost' filesep 'prePlacePost' filesep sourceFile.name]);
        end
    end
    cells_prePlace{clust_i,:} = setdiff(intersect(pp_cells,ppre_cells),cells_prePlacePost{clust_i,:});
    for cell_i = 1:length(cells_prePlace{clust_i,:})
        sourceFile = dir(['*_ROI_' num2str(cells_prePlace{clust_i}(cell_i)) '_cluster_' num2str(clust_i+1) '.tif']);
        if ~isempty(sourceFile)
            movefile(sourceFile.name,[pwd filesep 'cells_prePlacePost' filesep 'prePlace' filesep sourceFile.name]);
        end
    end
    cells_prePost{clust_i,:} = setdiff(intersect(ppre_cells,ppost_cells),cells_prePlacePost{clust_i,:});
    for cell_i = 1:length(cells_prePost{clust_i,:})
        sourceFile = dir(['*_ROI_' num2str(cells_prePost{clust_i}(cell_i)) '_cluster_' num2str(clust_i+1) '.tif']);
        if ~isempty(sourceFile)
            movefile(sourceFile.name,[pwd filesep 'cells_prePlacePost' filesep 'prePost' filesep sourceFile.name]);
        end
    end
    cells_placePost{clust_i,:} = setdiff(intersect(pp_cells,ppost_cells),cells_prePlacePost{clust_i,:});
    for cell_i = 1:length(cells_placePost{clust_i,:})
        sourceFile = dir(['*_ROI_' num2str(cells_placePost{clust_i}(cell_i)) '_cluster_' num2str(clust_i+1) '.tif']);
        if ~isempty(sourceFile)
            movefile(sourceFile.name,[pwd filesep 'cells_prePlacePost' filesep 'placePost' filesep sourceFile.name]);
        end
    end    
    cells_pre{clust_i,:} = setdiff(setdiff(ppre_cells,cells_prePlace{clust_i,:}),cells_prePost{clust_i,:});
    for cell_i = 1:length(cells_pre{clust_i,:})
        sourceFile = dir(['*_ROI_' num2str(cells_pre{clust_i}(cell_i)) '_cluster_' num2str(clust_i+1) '.tif']);
        if ~isempty(sourceFile)
            movefile(sourceFile.name,[pwd filesep 'cells_prePlacePost' filesep 'pre' filesep sourceFile.name]);
        end
    end
    cells_place{clust_i,:} = setdiff(setdiff(pp_cells,cells_placePost{clust_i,:}),cells_prePlace{clust_i,:});
    for cell_i = 1:length(cells_place{clust_i,:})
        sourceFile = dir(['*_ROI_' num2str(cells_place{clust_i}(cell_i)) '_cluster_' num2str(clust_i+1) '.tif']);
        if ~isempty(sourceFile)
            movefile(sourceFile.name,[pwd filesep 'cells_prePlacePost' filesep 'place' filesep sourceFile.name]);
        end
    end
    cells_post{clust_i,:} = setdiff(setdiff(ppost_cells,cells_placePost{clust_i,:}),cells_prePost{clust_i,:});
    for cell_i = 1:length(cells_post{clust_i,:})
        sourceFile = dir(['*_ROI_' num2str(cells_post{clust_i}(cell_i)) '_cluster_' num2str(clust_i+1) '.tif']);
        if ~isempty(sourceFile)
            movefile(sourceFile.name,[pwd filesep 'cells_prePlacePost' filesep 'post' filesep sourceFile.name]);
        end
    end
    
    
end

