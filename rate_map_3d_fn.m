function [D,place_field_density_raw_with_NaN_FE_smoothed,max_FR_smoothed_map,information_per_spike_3_D_FE_smoothed,sparsity_3_D_FE_smoothed,mat_timespent_density_raw]...
    = rate_map_3d_fn(pos,firpos,figure_flag,view_para)
%% Parameters for rate map
frames_per_second = 180;
room_dimensions = [0 5400 0 5400 0 2000]; % x,y,z
X_room_dimensions = room_dimensions(1:2);
Y_room_dimensions = room_dimensions(3:4);
Z_room_dimensions = room_dimensions(5:6);
bin_size_pixels = 100;
mult_factor = 5;
min_distance_from_trajectory = 300; % in mm: This is the minimal distance we allow a voxel to be from
% the an "Active" voxel, i.e., a voxel which had some video time in it.
Min_time_in_voxel = 1; % in sec
smoothing_voxel_size_cm = 50;
bin_size_pixels_2_D = 200;
h_rate_map = 1.5; % This is what Hafting2005 defines as the smoothing factor which is actually
%the std of the gaussian kernel.
sigmaa = h_rate_map; % This the standard deviation of the kernel, which determines its width.
% for larger values, the kernel and is in fact wider. Our 1.5 bins
% kernel is relatively conservative in comparison to what people
% usually use.
hsize =4*round(sigmaa)+1; % hsize is the size of the kernel: I define it as 5*round(sigma)+1
% Just in order to make sure it is long enough to facilitate suffficent
% fading of the gaussian. In practise, it is much longer then what we
% actually need!
gaussian_kernel = fspecial('gaussian',hsize,sigmaa);
time_spent_minimum = 0.2 ; % Will discard pixels in which the animal spent < this minimal time (seconds)
adaptive_smooth = 1;
%% Main loop for 3D rate map computation
x_video = pos(:,1);
y_video = pos(:,2);
z_video = pos(:,3);
x_spikes = firpos(:,1);
y_spikes = firpos(:,2);
z_spikes = firpos(:,3);
% x_spikes = x_spikes + 2600; y_spikes = y_spikes + 2600;
% x_video = x_video + 2600; y_video = y_video + 2600;

tic
mat_spike_density_raw = zeros( diff(Z_room_dimensions)/bin_size_pixels,...
    diff(Y_room_dimensions)/bin_size_pixels,diff(X_room_dimensions)/bin_size_pixels ) + NaN ; % Initialize
mat_timespent_density_raw =  zeros( diff(Z_room_dimensions)/bin_size_pixels,...
    diff(Y_room_dimensions)/bin_size_pixels,diff(X_room_dimensions)/bin_size_pixels )  ; % Initialize

place_field_density_raw =  zeros( diff(Z_room_dimensions)/bin_size_pixels,...
    diff(Y_room_dimensions)/bin_size_pixels,diff(X_room_dimensions)/bin_size_pixels ) + NaN ; % Initialize

Spike_IXs_per_voxel = cell( diff(Z_room_dimensions)/bin_size_pixels,...
    diff(Y_room_dimensions)/bin_size_pixels,diff(X_room_dimensions)/bin_size_pixels ) ; % Initialize

Video_IXs_per_voxel = cell( diff(Z_room_dimensions)/bin_size_pixels,...
    diff(Y_room_dimensions)/bin_size_pixels,diff(X_room_dimensions)/bin_size_pixels ) ; % Initialize

for ii_x_bin = 1 : diff(X_room_dimensions)/bin_size_pixels , % Loop over x-bins
    for ii_y_bin = 1 : diff(Y_room_dimensions)/bin_size_pixels , % Loop over y-bins
        for ii_z_bin = 1 : diff(Z_room_dimensions)/bin_size_pixels , % Loop over z-bins
            % Spike Density:
            mat_spike_density_raw( ii_z_bin, ii_y_bin,ii_x_bin) = ... % All the data
                sum( x_spikes >= 1 + bin_size_pixels*(ii_x_bin-1) & ...
                x_spikes <  1 + bin_size_pixels*(ii_x_bin) & ...
                y_spikes >= 1 + bin_size_pixels*(ii_y_bin-1) & ...
                y_spikes <  1 + bin_size_pixels*(ii_y_bin) & ...
                z_spikes >= 1 + bin_size_pixels*(ii_z_bin-1) & ...
                z_spikes <  1 + bin_size_pixels*(ii_z_bin)) ;
            
            Spike_IXs_per_voxel{ii_z_bin, ii_y_bin, ii_x_bin} = find( x_spikes >= 1 + bin_size_pixels*(ii_x_bin-1) & ...
                x_spikes <  1 + bin_size_pixels*(ii_x_bin) & ...
                y_spikes >= 1 + bin_size_pixels*(ii_y_bin-1) & ...
                y_spikes <  1 + bin_size_pixels*(ii_y_bin) & ...
                z_spikes >= 1 + bin_size_pixels*(ii_z_bin-1) & ...
                z_spikes<  1 + bin_size_pixels*(ii_z_bin)) ;
            
            
            % Time-Spent Density:
            mat_timespent_density_raw( ii_z_bin, ii_y_bin, ii_x_bin ) = ...
                sum( x_video >= 1 + bin_size_pixels*(ii_x_bin-1) & ...
                x_video <  1 + bin_size_pixels*(ii_x_bin) & ...
                y_video >= 1 + bin_size_pixels*(ii_y_bin-1) & ...
                y_video <  1 + bin_size_pixels*(ii_y_bin) & ...
                z_video >= 1 + bin_size_pixels*(ii_z_bin-1) & ...
                z_video <  1 + bin_size_pixels*(ii_z_bin) ) ;
            
            Video_IXs_per_voxel{ii_z_bin, ii_y_bin, ii_x_bin} = find( x_video >= 1 + bin_size_pixels*(ii_x_bin-1) & ...
                x_video <  1 + bin_size_pixels*(ii_x_bin) & ...
                y_video >= 1 + bin_size_pixels*(ii_y_bin-1) & ...
                y_video <  1 + bin_size_pixels*(ii_y_bin) & ...
                z_video >= 1 + bin_size_pixels*(ii_z_bin-1) & ...
                z_video <  1 + bin_size_pixels*(ii_z_bin) ) ;
            
            % Normalize Time-Spent Density from Video-Frames-Spent to Seconds-Spent
            
            mat_timespent_density_raw( ii_z_bin, ii_y_bin, ii_x_bin ) = ...
                mat_timespent_density_raw( ii_z_bin, ii_y_bin, ii_x_bin ) / frames_per_second ;
            
        end
    end
end

if ~adaptive_smooth
    place_field_density_raw_with_NaN = mat_spike_density_raw./mat_timespent_density_raw;
    mat_time_density_raw_with_NaN = mat_timespent_density_raw;
else
    
    place_field_density_raw_with_NaN = zeros( diff(Z_room_dimensions)/bin_size_pixels,...
        diff(Y_room_dimensions)/bin_size_pixels,diff(X_room_dimensions)/bin_size_pixels ) + NaN ; % Initialize
    mat_time_density_raw_with_NaN = zeros( diff(Z_room_dimensions)/bin_size_pixels,...
        diff(Y_room_dimensions)/bin_size_pixels,diff(X_room_dimensions)/bin_size_pixels ) + NaN ; % Initialize
    mat_spike_count_raw_with_NaN = zeros( diff(Z_room_dimensions)/bin_size_pixels,...
        diff(Y_room_dimensions)/bin_size_pixels,diff(X_room_dimensions)/bin_size_pixels ) + NaN ; % Initialize
    
    % Create a MUCH larger "fictive" room into which our real room will go
    % into (this will be used later for the summation of the timespent).
    mat_timespent_density_raw_fictive = zeros( mult_factor*diff(Z_room_dimensions)/bin_size_pixels,...
        mult_factor*diff(Y_room_dimensions)/bin_size_pixels,mult_factor*diff(X_room_dimensions)/bin_size_pixels )  ; % Initialize
    mat_spike_density_raw_fictive = zeros( mult_factor*diff(Z_room_dimensions)/bin_size_pixels,...
        mult_factor*diff(Y_room_dimensions)/bin_size_pixels,mult_factor*diff(X_room_dimensions)/bin_size_pixels )  ; % Initialize
    
    % Note below the matrixs are of the REAL size as we will be
    % inserting the data directly into them:
    mat_associated_video_IXs = cell( diff(Z_room_dimensions)/bin_size_pixels,...
        diff(Y_room_dimensions)/bin_size_pixels,diff(X_room_dimensions)/bin_size_pixels ) ; % Initialize
    mat_associated_spike_IXs = cell( diff(Z_room_dimensions)/bin_size_pixels,...
        diff(Y_room_dimensions)/bin_size_pixels,diff(X_room_dimensions)/bin_size_pixels ) ; % Initialize
    
    
    real_Z_IXs = (mult_factor/2 - 0.5)*diff(Z_room_dimensions)/bin_size_pixels + 1:(mult_factor/2 + 0.5)*diff(Z_room_dimensions)/bin_size_pixels;
    real_Y_IXs = (mult_factor/2 - 0.5)*diff(Y_room_dimensions)/bin_size_pixels + 1:(mult_factor/2 + 0.5)*diff(Y_room_dimensions)/bin_size_pixels;
    real_X_IXs = (mult_factor/2 - 0.5)*diff(X_room_dimensions)/bin_size_pixels + 1:(mult_factor/2 + 0.5)*diff(X_room_dimensions)/bin_size_pixels;
    
    % insert our real maxtrix into the fictive one:
    mat_timespent_density_raw_fictive(real_Z_IXs,real_Y_IXs,real_X_IXs) = mat_timespent_density_raw;
    mat_spike_density_raw_fictive(real_Z_IXs,real_Y_IXs,real_X_IXs) = mat_spike_density_raw;
    
    % Now for each voxel we would like to compute the rate by using
    % adaptive smoothing implimented by increase the size of the voxel
    % until we have over a certain number of timespent in it (defined via
    % the variable: Min_time_in_voxel)
    validation_neighbors = min_distance_from_trajectory/bin_size_pixels;
    for ii_x_bin = 1 : diff(X_room_dimensions)/bin_size_pixels  , % Loop over x-bins
        for ii_y_bin = 1 : diff(Y_room_dimensions)/bin_size_pixels , % Loop over y-bins
            for ii_z_bin = 1 : diff(Z_room_dimensions)/bin_size_pixels , % Loop over z-bins
                current_voxel_timespent = mat_timespent_density_raw_fictive(real_Z_IXs(ii_z_bin),real_Y_IXs(ii_y_bin),real_X_IXs(ii_x_bin));
                neighbors_counter = 0;
                while (current_voxel_timespent<=Min_time_in_voxel)
                    neighbors_counter = neighbors_counter + 1; % Increase the voxel by one in each dimension
                    current_time_voxel = mat_timespent_density_raw_fictive(real_Z_IXs(ii_z_bin)-neighbors_counter:real_Z_IXs(ii_z_bin)+neighbors_counter,...
                        real_Y_IXs(ii_y_bin)-neighbors_counter:real_Y_IXs(ii_y_bin)+neighbors_counter,...
                        real_X_IXs(ii_x_bin)-neighbors_counter:real_X_IXs(ii_x_bin)+neighbors_counter);
                    current_voxel_timespent = sum(sum(sum(current_time_voxel)));
                end
                
                current_voxel_IX_relative_to_original_mat =[ii_z_bin-neighbors_counter:ii_z_bin+neighbors_counter;...
                    ii_y_bin-neighbors_counter:ii_y_bin+neighbors_counter;...
                    ii_x_bin-neighbors_counter:ii_x_bin+neighbors_counter];
                
                
                % Check if this voxel is anywhere near an "active"
                % voxel
                current_time_valid = mat_timespent_density_raw_fictive(real_Z_IXs(ii_z_bin)-validation_neighbors:real_Z_IXs(ii_z_bin)+validation_neighbors,...
                    real_Y_IXs(ii_y_bin)-validation_neighbors:real_Y_IXs(ii_y_bin)+validation_neighbors,...
                    real_X_IXs(ii_x_bin)-validation_neighbors:real_X_IXs(ii_x_bin)+validation_neighbors);
                
                %                      current_time_valid_DBG = mat_timespent_density_raw(ii_z_bin-validation_neighbors:ii_z_bin,...
                %                         ii_y_bin-validation_neighbors:ii_y_bin,...
                %                         ii_x_bin-validation_neighbors:ii_x_bin);
                %
                %                     current_time_valid_DBG = mat_timespent_density_raw(ii_z_bin-validation_neighbors,...
                %                         ii_y_bin-validation_neighbors,...
                %                         ii_x_bin);
                
                current_spike_voxel = mat_spike_density_raw_fictive(real_Z_IXs(ii_z_bin)-neighbors_counter:real_Z_IXs(ii_z_bin)+neighbors_counter,...
                    real_Y_IXs(ii_y_bin)-neighbors_counter:real_Y_IXs(ii_y_bin)+neighbors_counter,...
                    real_X_IXs(ii_x_bin)-neighbors_counter:real_X_IXs(ii_x_bin)+neighbors_counter);
                
                current_voxel_spike_count = sum(sum(sum(current_spike_voxel)));
                
                if (sum(sum(sum(current_time_valid)))== 0) % i.e, there was NO activity in this distance from the central voxel (including it of course)
                    current_voxel_timespent = NaN; % we will not be using this voxel for analysis give it a vlaue of NaN here.
                    current_spike_voxel = NaN;
                else
                    % save the relevant position of animal + spike relevant
                    % for the computation in this voxel
                    for ii_x = current_voxel_IX_relative_to_original_mat(3,:)
                        for ii_y = current_voxel_IX_relative_to_original_mat(2,:)
                            for ii_z = current_voxel_IX_relative_to_original_mat(1,:)
                                current_IX_position = [ii_z,ii_y,ii_x];
                                % Check if Z IX is range
                                check_z = ((ii_z>0)&(ii_z<size(mat_associated_video_IXs,1)));
                                check_y = ((ii_y>0)&(ii_y<size(mat_associated_video_IXs,2)));
                                check_x = ((ii_x>0)&(ii_x<size(mat_associated_video_IXs,3)));
                                if (sum([check_z,check_y,check_x]) == 3)% meaning all are positive IXs
                                    mat_associated_video_IXs{ii_z_bin,ii_y_bin,ii_x_bin} = ...
                                        [mat_associated_video_IXs{ii_z_bin,ii_y_bin,ii_x_bin},...
                                        (Video_IXs_per_voxel{current_IX_position(1),...
                                        current_IX_position(2),current_IX_position(3)})'];
                                    
                                    mat_associated_spike_IXs{ii_z_bin,ii_y_bin,ii_x_bin} = ...
                                        [mat_associated_spike_IXs{ii_z_bin,ii_y_bin,ii_x_bin},...
                                        (Spike_IXs_per_voxel{current_IX_position(1),...
                                        current_IX_position(2),current_IX_position(3)})'];
                                else end
                                
                            end
                        end
                    end
                end
                place_field_density_raw_with_NaN(ii_z_bin,...
                    ii_y_bin,...
                    ii_x_bin) = current_voxel_spike_count/current_voxel_timespent;
                mat_time_density_raw_with_NaN(ii_z_bin,...
                    ii_y_bin,...
                    ii_x_bin) = current_voxel_timespent;
                mat_spike_count_raw_with_NaN(ii_z_bin,...
                    ii_y_bin,...
                    ii_x_bin)= current_voxel_spike_count;
                
            end
        end
    end
end
%% Plot the 3D rate map
% place_field_density_raw_with_NaN_reshape = permute(place_field_density_raw_with_NaN,[2,3,1]);
% figure
% vol3d('cdata',place_field_density_raw_with_NaN_reshape,'texture','3D');
% % h = slice(place_field_density_raw_with_NaN, [], [], 1:size(place_field_density_raw_with_NaN,3));
% % set(h, 'EdgeColor','none', 'FaceColor','interp')
% % alpha(.05)
% % colormap(parula)
% view([-81 10])
%
transperancy_thres = 20; % in percentages of maximal firing rate
transperancy_val_of_mx_FR = 4; % in percentages of maximal firing rate
% along each dimension. The reason why we do this is to exclude cases where
% a single voxel outlier would affect the field size estimation.
X_room_dimensions_cm = room_dimensions(1:2)/10;
Y_room_dimensions_cm = room_dimensions(3:4)/10;
Z_room_dimensions_cm = room_dimensions(5:6)/10;

color_list = [ 1 0 0 ; ... % List of 18 colors for plotting spikes of different units
    0 0.7 0 ; ...
    0 0 1 ; ...
    0.5 0.0 0.5 ; ...
    0.5 0.5 0.0 ; ...
    0 0.5 0.5 ; ...
    0.9 0.1 0.5 ; ...
    0.9 0.5 0.1 ; ...
    0.1 0.9 0.5 ; ...
    0.1 0.5 0.9 ; ...
    0.5 0.9 0.1 ; ...
    0.5 0.1 0.9 ; ...
    0.3 0.7 0.7 ; ...
    0.7 0.3 0.7 ; ...
    0.7 0.7 0.3 ; ...
    0.3 0.3 0.7 ; ...
    0.3 0.7 0.3 ; ...
    0.7 0.3 0.3 ] ;

max_FR = max(max(max(place_field_density_raw_with_NaN)));
transperancy_thres_val = transperancy_thres*max_FR/100;
value_set_at_below_thres_voxels = transperancy_val_of_mx_FR*max_FR/100;
IXs = find(place_field_density_raw_with_NaN<transperancy_thres_val);
place_field_density_raw_with_NaN_FE_adjusted = place_field_density_raw_with_NaN;
if ~isempty(IXs)
    place_field_density_raw_with_NaN_FE_adjusted(IXs) = value_set_at_below_thres_voxels;
    place_field_density_raw_with_NaN_FE_adjusted(IXs(1))=0; % To set the scale from zero
else end
place_field_density_raw_with_NaN_FE_adjusted_reshaped = permute(place_field_density_raw_with_NaN_FE_adjusted,[2,3,1]);
%D = squeeze(place_field_density_raw_with_NaN_FE_adjusted_reshaped);
D = (place_field_density_raw_with_NaN_FE_adjusted_reshaped);

idx_NaN_temp = find(isnan(place_field_density_raw_with_NaN));
place_field_density_raw_FE_smoothed_no_NaN = place_field_density_raw_with_NaN;
place_field_density_raw_FE_smoothed_no_NaN(idx_NaN_temp) = 0;
mat_time_density_raw_FE_smoothed_no_NaN = mat_time_density_raw_with_NaN;
mat_time_density_raw_FE_smoothed_no_NaN(idx_NaN_temp) = 0;
% if rem(size(place_field_density_raw_FE_smoothed_no_NaN,1),2) || rem(size(place_field_density_raw_FE_smoothed_no_NaN,2),2) ...
%         || rem(size(place_field_density_raw_FE_smoothed_no_NaN,3),2) == 0
%     place_field_density_raw_FE_smoothed_no_NaN = place_field_density_raw_FE_smoothed_no_NaN(1:end-1,1:end-1,1:end-1);
%     mat_time_density_raw_FE_smoothed_no_NaN = mat_time_density_raw_FE_smoothed_no_NaN(1:end-1,1:end-1,1:end-1);
% end

% We want to smooth the 3D matrix WITHOUT the NaNs:
place_field_density_raw_with_NaN_FE_smoothed = smooth3(place_field_density_raw_FE_smoothed_no_NaN,...
    'gaussian',[smoothing_voxel_size_cm/10,smoothing_voxel_size_cm/10,smoothing_voxel_size_cm/10]);
mat_time_density_raw_with_NaN_FE_smoothed = smooth3(mat_time_density_raw_FE_smoothed_no_NaN,...
    'gaussian',[smoothing_voxel_size_cm/10,smoothing_voxel_size_cm/10,smoothing_voxel_size_cm/10]);

% Now return the NaNs so that we won't count the unvisited voxels:
place_field_density_raw_with_NaN_FE_smoothed(idx_NaN_temp) = NaN;
mat_time_density_raw_with_NaN_FE_smoothed(idx_NaN_temp) = NaN;


% -----  INFORMATION PER SPIKE  (computed for the SMOOTHED field): -----
% Information_per_spike = sum( p_i * ( r_i / r ) * log2( r_i / r ) )
%    Where:
%       r_i = firing rate in bin i ;
%       p_i = occupancy of bin i = time-spent by bat in bin i / total time spent in all bins ;
%       r = mean( r_i ) = overall mean firing rate (mean over all the pixels)
% See: Skaggs WE, McNaughton BL, Wilson MA, Barnes CA, Hippocampus 6(2), 149-172 (1996).

%find the not-NaN IXs;
idx_notNaN_PlaceField_FE = find(~isnan(place_field_density_raw_with_NaN_FE_smoothed));
% For the smoothed map:
r_i = place_field_density_raw_with_NaN_FE_smoothed( idx_notNaN_PlaceField_FE ); % Use the SMOOTHED Place Field
p_i = mat_time_density_raw_with_NaN_FE_smoothed( idx_notNaN_PlaceField_FE ) ./ ...
    sum( mat_time_density_raw_with_NaN_FE_smoothed( idx_notNaN_PlaceField_FE ) ) ;
r_i = r_i(:) ; p_i = p_i(:) ; % Turn r_i and p_i into Column vectors
% % %             r = mean( r_i ) ;
r = sum( r_i .* p_i );
information_per_spike_3_D_FE_smoothed = sum( p_i .* ( r_i / r ) .* log2( ( r_i + eps ) / r ) ) ; % I added a tiny number to avoid log(0)

% -----  SPARSITY  (computed for the SMOOTHED field): -----
% Sparsity = <r_i>^2 / <r_i^2> = sum( p_i * r_i )^2 / sum( p_i * r_i^2 )
%    Where:
%       r_i = firing rate in bin i ;
%       p_i = occupancy of bin i = time-spent by bat in bin i / total time spent in all bins.
% See: Skaggs WE, McNaughton BL, Wilson MA, Barnes CA, Hippocampus 6(2), 149-172 (1996).
sparsity_3_D_FE_smoothed = sum( p_i .* r_i )^2 / sum( p_i .* ( r_i .^2 ) ) ; % Sparsity

max_FR = max(max(max(place_field_density_raw_with_NaN_FE_smoothed)));
transperancy_thres_val = transperancy_thres*max_FR/100;
value_set_at_below_thres_voxels = transperancy_val_of_mx_FR*max_FR/100;
IXs = find(place_field_density_raw_with_NaN_FE_smoothed<transperancy_thres_val);
place_field_density_raw_with_NaN_FE_smoothed_adjusted = place_field_density_raw_with_NaN_FE_smoothed;
if ~isempty(IXs)
    place_field_density_raw_with_NaN_FE_smoothed_adjusted(IXs) = value_set_at_below_thres_voxels;
    place_field_density_raw_with_NaN_FE_smoothed_adjusted(IXs(1))=0; % To set the scale from zero
else
end
% ax_smoothed_PF = axes('position',[.75  .7  .17  .19]);
colormap_map = jet ;
colormap_map(64,:) = [1 1 1] ; % Set the highest-value pixel to White = [1 1 1] : I will use White as the Color of NaN pixels
cdata_mat = place_field_density_raw_with_NaN_FE_smoothed / max(max(max(place_field_density_raw_with_NaN_FE_smoothed))) * ...
    ( size(colormap_map,1) - 1 ) ; % Scale the colors in 'cdata' to ( 1 : length(colormap) - 1 ), so that the highest value will be reserved to NaN's
NaN_IX = find(isnan( place_field_density_raw_with_NaN_FE_smoothed));
cdata_mat( NaN_IX ) = size(colormap_map,1) ; % Replace NaN values with the highest values (they will be colored white)
place_field_density_raw_with_NaN_FE_smoothed_reshaped = permute(place_field_density_raw_with_NaN_FE_smoothed_adjusted,[2,3,1]);
D = squeeze(place_field_density_raw_with_NaN_FE_smoothed_reshaped);
max_FR_smoothed_map = max(max(max(place_field_density_raw_with_NaN_FE_smoothed)));

%Plot the 3D rate map
if figure_flag == 1
    figure('pos',[550 300 1300 800])
    pos = [x_video y_video z_video];
    pos_norm(:,1) = (pos(:,1)-min(pos(:,1)))/(max(pos(:,1))-min(pos(:,1)));
    pos_norm(:,2) = (pos(:,2)-min(pos(:,2)))/(max(pos(:,2))-min(pos(:,2)));
    pos_norm(:,3) = (pos(:,3)-min(pos(:,3)))/(max(pos(:,3))-min(pos(:,3)));
    b = D;
    pos_rescale(:,1) = (size(b,1)-1)*pos_norm(:,1) +1;
    pos_rescale(:,2) = (size(b,2)-1)*pos_norm(:,2) +1;
    pos_rescale(:,3) = (size(b,3)-1)*pos_norm(:,3) +1;
    % plot3(pos_rescale(:,1),pos_rescale(:,2),pos_rescale(:,3),'color', [0.5 0.5 0.5],'Linewidth',0.5);
    
    cuboid_dim = [max(pos_rescale(:,1))+2 max(pos_rescale(:,2)) max(pos_rescale(:,3))+5]; % in meter
    CV = [mean(pos_rescale(:,1));mean(pos_rescale(:,2));mean(pos_rescale(:,3))];
    EA = [0;0;0];   colr = [0.5 0.5 0.5];   alph = 0.01;
    [CuboidHandle, verts, facs] = DrawCuboid(cuboid_dim,CV,EA,colr,alph); hold on
    
    
    h = vol3d('cdata',D,'texture','2D');
    view(3);
    % Update view since 'texture' = '2D'
    vol3d(h);
    axis tight; axis equal;
    axis off
    %daspect([1 1 .4])
    alphamap('rampup');
    alphamap(.6 .* alphamap);
    grid on
    set(gca,'XTick',[0,20,40,60],'XTickLabel',([0,200,400,600]))
    set(gca,'YTick',[0,25,50],'YTickLabel',([0,250,500]))
    set(gca,'ZTick',[0,10,20,30],'ZTickLabel',([0,100,200,300]))
    xlabel('X(cm)')
    ylabel('Y(cm)')
    zlabel('Z(cm)')
    
    % title({['3-D Rate Map - Smoothed'],['3-D Gaussian kernel size = ', num2str(smoothing_voxel_size_cm),' cm'],...
    %     ['Information per spike - smoothed rate map = ', num2str(information_per_spike_3_D_FE_smoothed,2)],...
    %     ['Sparsity = ', num2str(sparsity_3_D_FE_smoothed,2)]})
    
    text(-5,1,32,...
        strcat(num2str(max_FR_smoothed_map,3),{' '},'Hz'),'FontSize',40,'FontWeight','bold')
    toc
    view([view_para(1) view_para(2)])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
    colormap(jet)
end
% axis([0 5200 0 5200 0 2200]) % room dimensions
% grid on
% set(gca,'XTick',[0,2000,4000,6000],'XTickLabel',([0,2000,4000,6000]),'xlim',[0,5200])
% set(gca,'YTick',[0,2500,5000],'YTickLabel',([0,2500,5000]),'ylim',[0,5200])
% set(gca,'ZTick',[0,1000,2000],'ZTickLabel',([0,1000,2000]),'zlim',[0,2000])
end