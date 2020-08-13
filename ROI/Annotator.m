classdef Annotator < handle
    %ANNOTATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % image
        name
        file
        
        % annotation file
        annot_file = '';
    end
    
    properties (Access=protected)
        % image
        image
        
        % dimensions
        width
        height
        
        % annotations
        annotations = [];
        
        % mode
        saved = true;
        
        % handles
        win
        axes
        
        % gui elements
        gui_toolbar
        
        % plot handles
        plot_annotations
        
        % cached masks
        cached_masks
    end
    
    methods
        function AN = Annotator(file, annot_file)
            % make full path
            file = fullfile(pwd, file);
            
            if ~exist(file, 'file')
                error('Invalid image file.');
            end
            disp(file);
            
            % extract parts
            [path, nm] = fileparts(file);
            
            % load image
            image = imread(file);
            
            % set parameters
            AN.name = nm;
            AN.file = file;
            AN.image = image;
            
            % store images
            AN.width = size(image, 2);
            AN.height = size(image, 1);
            
            % make color
            color_bg = [0.85 0.85 0.85];
            
            % get screen size
            screen = get(0, 'ScreenSize');
            
            % inital dimensions
            h = size(image, 1) ;
            w = size(image, 2);
            x = max((screen(3) - w) / 2, 0);
            y = max((screen(2) - h) / 2, 0);
            
            % create viewer window
            AN.win = figure('Visible', 'on', 'Name', nm, ...
                'Position', [x y w h], 'NumberTitle', 'off', 'Toolbar', ...
                'none', 'MenuBar', 'none', 'Resize', 'off', 'Color', ...
                color_bg);
            
            % set 
            set(AN.win, 'PaperPositionMode', 'auto');
            set(AN.win, 'InvertHardcopy', 'off');
            set(AN.win, 'Units', 'pixels');
            set(AN.win, 'CloseRequestFcn', @AN.cb_beforeCloseWindow);
            set(AN.win, 'DeleteFcn', @AN.cb_closeWindow);
            
            % toolbar
            AN.gui_toolbar = uitoolbar('Parent', AN.win);
            
            % add new button
            [ico, ~, alpha] = imread(fullfile(matlabroot, 'toolbox', 'matlab', 'icons', 'file_new.png'));
            if isa(ico, 'uint8')
                ico = double(ico) / (256 - 1);
            elseif isa(ico, 'uint16')
                ico = double(ico) / (256 * 256 - 1);
            end
            ico(repmat(alpha == 0, 1, 1, size(ico, 3))) = nan;
            uipushtool('Parent', AN.gui_toolbar, 'CData', ico, ...
                'ClickedCallback', @AN.cb_new, 'TooltipString', ...
                'New');
            
            % add open button
            [ico, ~, alpha] = imread(fullfile(matlabroot, 'toolbox', 'matlab','icons', 'file_open.png'));
            if isa(ico, 'uint8')
                ico = double(ico) / (256 - 1);
            elseif isa(ico, 'uint16')
                ico = double(ico) / (256 * 256 - 1);
            end
            ico(repmat(alpha == 0, 1, 1, size(ico, 3))) = nan;
            uipushtool('Parent', AN.gui_toolbar, 'CData', ico, ...
                'ClickedCallback', @AN.cb_load, 'TooltipString', ...
                'Open');
            
            % add save button
            [ico, ~, alpha] = imread(fullfile(matlabroot, 'toolbox', 'matlab', 'icons', 'file_save.png'));
            if isa(ico, 'uint8')
                ico = double(ico) / (256 - 1);
            elseif isa(ico, 'uint16')
                ico = double(ico) / (256 * 256 - 1);
            end
            ico(repmat(alpha == 0, 1, 1, size(ico, 3))) = nan;
            uipushtool('Parent', AN.gui_toolbar, 'CData', ico, ...
                'ClickedCallback', @AN.cb_save, 'TooltipString', ...
                'Save');
            
            % get axes
            AN.axes = axes('Parent', AN.win);
            axis off;
            
            % show image
            imshow(AN.image, 'Parent', AN.axes, 'Border', 'tight');
            pan off;
            
            % auto load annotations if specified or same file name exists
            if exist('annot_file', 'var') && ~isempty(annot_file)
                AN.loadAnnotations(annot_file);
            else
                default_annot_file = [path filesep nm '.mat'];
                if exist(default_annot_file, 'file')
                    AN.loadAnnotations(default_annot_file);
                end
            end
        end
        
        function delete(AN)
            try
                delete(AN.win);
            catch err %#ok<NASGU>
            end
        end
        
        function cb_new(AN, h, event)
            if length(AN.plot_annotations) ~= size(AN.annotations, 1)
                error('Annotion list out of sync. Try reopening annotator.');
            end
            
            % get annotation
            h = AN.drawAnnotation();
            if ~isvalid(h)
                return;
            end
            
            % mark as unsaved
            AN.saved = false;
            
            % add to list
            AN.plot_annotations{end + 1} = h;
            AN.annotations = [AN.annotations; h.getPosition()]; 
        end
        
        function cb_load(AN, h, event)
            [filename, pathname] = uigetfile({'*.mat', 'MATLAB File (*.mat)'; '*.*', 'All Files'}, 'Load annotations');
            
            % was anceled?
            if isequal(filename, 0) || isequal(pathname, 0)
                return;
            end
            
            % load file
            AN.loadAnnotations(fullfile(pathname, filename));
        end
        
        function cb_save(AN, h, event)
            % figure out default name
            if isempty(AN.annot_file)
                [path, nm, ~] = fileparts(AN.file);
                def_name = [path filesep nm '.mat'];
            else
                def_name = AN.annot_file;
            end
            
            % show save window
            [filename, pathname] = uiputfile({'*.mat', 'MATLAB File (*.mat)'; '*.*', 'All Files'}, 'Save annotations', def_name);
            
            % was anceled?
            if isequal(filename, 0) || isequal(pathname, 0)
                return;
            end
            
            % save
            AN.saveAnnotations(fullfile(pathname, filename));
        end
        
        function cb_beforeCloseWindow(AN, h, event)
            % clean annotations
            AN.cleanAnnotations();
            
            % cache masks
            AN.cached_masks = AN.getMasks();
            
            % is unsaved?
            if ~AN.saved
                % prompt to save
                response = questdlg('Do you want to save changes before closing?', 'Save Changes', 'Cancel', 'No', 'Yes', 'Yes');
                switch response
                    case 'Yes'
                        AN.cb_save(h, event);
                    case 'Cancel'
                        return
                end
            end
            
            % clear image
            delete(gcf);
        end
        
        function cb_closeWindow(AN, h, event)
            % clear image
            clear AN.image;
        end
        
        function loadAnnotations(AN, fl)
            % load file
            d = load(fl);
            
            % check file
            if ~isfield(d, 'annotations')
                warning('Invalid annotations file.');
                return
            end
            
            if d.width ~= AN.width && d.height ~= AN.height
                error('The annotation file has different dimensions.');
            end
            
            if ~strcmp(d.file, AN.file)
                warning('Annotations were potentially for a different image.');
            end
            
            % store file name
            AN.annot_file = fl;
            
            % copy data
            AN.annotations = d.annotations;
            
            % redraw
            AN.redrawAnnotations();
            
            % mark saved
            AN.saved = true;
        end
        
        function saveAnnotations(AN, fl)
            % clean annotations
            AN.cleanAnnotations();
            
            % extract variables
            name = AN.name; %#ok<NASGU,PROPLC>
            file = AN.file; %#ok<NASGU,PROPLC>
            %image = AN.image; %#ok<NASGU,PROPLC>
            annotations = AN.annotations; %#ok<NASGU,PROPLC>
            width = AN.width; %#ok<NASGU,PROPLC>
            height = AN.height; %#ok<NASGU,PROPLC>
            %masks = AN.cached_masks; %#ok<NASGU,PROPLC>
            % do save
            save(fl, '-v7.3', 'name', 'file', 'annotations', 'width', 'height');
            
            % store file name
            AN.annot_file = fl;
            
            % mark saved
            AN.saved = true;
        end
        
        function results = wait(AN)
            % block
            uiwait(AN.win);
            results = AN.annotations;
        end
        
        function results = getAnnotations(AN)
            AN.cleanAnnotations();
            results = AN.annotations;
        end
        
        function masks = getMasks(AN)
            if ~isvalid(AN.win)
                masks = AN.cached_masks;
                return;
            end
            
            AN.cleanAnnotations();
            
            ret = {};
            for i = 1:length(AN.plot_annotations)
                h = AN.plot_annotations{i};
                if ~isvalid(h)
                    continue;
                end
                
                ret{end + 1} = h.createMask();
            end
            masks = cat(3, ret{:});
        end
    end
    
    methods (Access=protected)
        function h = drawAnnotation(AN, position)
            if ~exist('position', 'var') || isempty(position)
                h = imellipse(AN.axes);
                if ~isvalid(h)
                    return;
                end
            else
                h = imellipse(AN.axes, position);
            end
            addNewPositionCallback(h, @(p) AN.annotationPositionUpdate(h, p));
        end
        
        function annotationPositionUpdate(AN, h, position)
            for i = 1:length(AN.plot_annotations)
                if AN.plot_annotations{i} == h
                    AN.annotations(i, :) = position;
                    AN.saved = false;
                    break;
                end
            end
        end
        
        function redrawAnnotations(AN)
            % hold axes
            hold(AN.axes, 'on');
            
            % remove existing plot
            if ~isempty(AN.plot_annotations)
                for i = 1:length(AN.plot_annotations)
                    h = AN.plot_annotations{i};
                    delete(h);
                end
                AN.plot_annotations = {};
            end
            
            % add new plot
            if ~isempty(AN.annotations)
                for i = 1:size(AN.annotations, 1)
                    % draw new annotation
                    h = AN.drawAnnotation(AN.annotations(i, :));
                    AN.plot_annotations{i} = h;
                end
            end
            
            % unhold axes
            hold(AN.axes, 'off');
        end
        
        function cleanAnnotations(AN)
            % only if window is still valid
            if ~isvalid(AN.win)
                return;
            end
            
            new_annotations = [];
            new_plot_annotations = {};
            for i = 1:length(AN.plot_annotations)
                % not valid
                if ~isvalid(AN.plot_annotations{i})
                    continue;
                end
                
                % add annotations
                new_annotations = [new_annotations; AN.annotations(i, :)];
                new_plot_annotations{end + 1} = AN.plot_annotations{i};
            end
            AN.annotations = new_annotations;
            AN.plot_annotations = new_plot_annotations;
        end
    end
end
