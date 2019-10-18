function [place_field_density_smoothed_with_NaN_X_Y,place_field_density_smoothed_with_NaN_Z_Y,place_field_density_smoothed_with_NaN_X_Z,...
    peak_firing_rate_X_Y,peak_firing_rate_Z_Y,peak_firing_rate_X_Z,SI_XY,SI_YZ,SI_XZ]...
    = rate_map_2d_fn(pos,firpos,figure_flag)
%% Parameters for rate map
frames_per_second = 120;
room_dimensions = [0 5400 0 5400 0 2000]; % x,y,z
% room_dimensions = [0 6000 0 6000 0 2000]; % x,y,z
X_room_dimensions = room_dimensions(1:2);
Y_room_dimensions = room_dimensions(3:4);
Z_room_dimensions = room_dimensions(5:6);
bin_size_pixels = 100;
mult_factor = 5;
min_distance_from_trajectory = 250; % in mm: This is the minimal distance we allow a voxel to be from
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
%% Main loop for 3D rate map computation
x_video = pos(:,1);
y_video = pos(:,2);
z_video = pos(:,3);
x_spikes = firpos(:,1);
y_spikes = firpos(:,2);
z_spikes = firpos(:,3);
% x_spikes = x_spikes + 2600; y_spikes = y_spikes + 2600;
% x_video = x_video + 2600; y_video = y_video + 2600;

%% Rate map XY
mat_spike_density_raw_X_Y = zeros( diff(Y_room_dimensions)/bin_size_pixels_2_D,...
    diff(X_room_dimensions)/bin_size_pixels_2_D ) + NaN ; % Initialize
mat_timespent_density_raw_X_Y =  zeros( diff(Y_room_dimensions)/bin_size_pixels_2_D,...
    diff(X_room_dimensions)/bin_size_pixels_2_D ) + NaN ; % Initialize
place_field_density_raw_X_Y =  zeros( diff(Y_room_dimensions)/bin_size_pixels_2_D,...
    diff(X_room_dimensions)/bin_size_pixels_2_D ) + NaN ; % Initialize

for ii_x_bin = 1 : diff(X_room_dimensions)/bin_size_pixels_2_D , % Loop over x-bins
    for ii_y_bin = 1 : diff(Y_room_dimensions)/bin_size_pixels_2_D , % Loop over y-bins
        % Spike Density:
        mat_spike_density_raw_X_Y( ii_y_bin,ii_x_bin) = ... % All the data
            sum( x_spikes >= 1 + bin_size_pixels_2_D*(ii_x_bin-1) & ...
            x_spikes <  1 + bin_size_pixels_2_D*(ii_x_bin) & ...
            y_spikes >= 1 + bin_size_pixels_2_D*(ii_y_bin-1) & ...
            y_spikes <  1 + bin_size_pixels_2_D*(ii_y_bin)) ;
        % Time-Spent Density:
        mat_timespent_density_raw_X_Y( ii_y_bin, ii_x_bin ) = ...
            sum( x_video >= 1 + bin_size_pixels_2_D*(ii_x_bin-1) & ...
            x_video <  1 + bin_size_pixels_2_D*(ii_x_bin) & ...
            y_video >= 1 + bin_size_pixels_2_D*(ii_y_bin-1) & ...
            y_video <  1 + bin_size_pixels_2_D*(ii_y_bin) ) ;
        % Normalize Time-Spent Density from Video-Frames-Spent to Seconds-Spent
        mat_timespent_density_raw_X_Y( ii_y_bin, ii_x_bin ) = ...
            mat_timespent_density_raw_X_Y( ii_y_bin, ii_x_bin ) / frames_per_second ;
        %             Discard pixels in which the animal Spent less than a certain Minimal amount of time --
        %             (this is computed for the "idx_include_VT" data only, usually resulting in
        % DIFFERENT pixels being discarded for the Full data):
        if ( mat_timespent_density_raw_X_Y( ii_y_bin, ii_x_bin) < time_spent_minimum ),
            mat_timespent_density_raw_X_Y( ii_y_bin, ii_x_bin) = 0 ; % Discard this time-spent-density pixel
            mat_timespent_density_raw_X_Y( ii_y_bin, ii_x_bin ) = 0 ; % Discard this spike-density pixel
            
        end
        if(mat_timespent_density_raw_X_Y~=0 & mat_timespent_density_raw_X_Y ==0)
            disp(2);
        end
    end
end

% Place Field = Spike Density / Time-Spent Density :
%                     warning off MATLAB:divideByZero ;
place_field_density_raw_X_Y = mat_spike_density_raw_X_Y ./ mat_timespent_density_raw_X_Y;
%                     warning on MATLAB:divideByZero ;


% Smoothing = convolve with gaussian kernel:
mat_spike_density_smoothed_X_Y = imfilter(mat_spike_density_raw_X_Y,gaussian_kernel);
mat_timespent_density_smoothed_X_Y = imfilter(mat_timespent_density_raw_X_Y,gaussian_kernel);


% Place Field smoothed = Spike Density smoothed / Time-Spent Density smoothed :
%                     warning off MATLAB:divideByZero ;
place_field_density_smoothed_X_Y = mat_spike_density_smoothed_X_Y ./ mat_timespent_density_smoothed_X_Y ;
%                     warning on MATLAB:divideByZero ;

% ======= Compute the PF density with NaN's at unvisited location (will later be presented as white bins in the PF figure) : ==========

% "Legalize" a bin (remove NaN) if the bat visited any of the bin's 8 closest neighbours:
%                     warning off all
idx_timespent_density = zeros( diff(Y_room_dimensions)/bin_size_pixels_2_D, diff(X_room_dimensions)/bin_size_pixels_2_D ) + NaN ; % Initialize
for ii_x_bin = 2 : (diff(X_room_dimensions)/bin_size_pixels_2_D - 1) , % Loop over x-bins, NOT INCL. CAMERA-VIEW EDGES
    for ii_y_bin = 2 : (diff(Y_room_dimensions)/bin_size_pixels_2_D - 1) , % Loop over y-bins, NOT INCL. CAMERA-VIEW EDGES
        matrix_3x3_of_neighbors = ...
            mat_timespent_density_raw_X_Y( ii_y_bin-1 : ii_y_bin+1, ii_x_bin-1 : ii_x_bin+1 ) ;
        sum_including_the_central_bin = sum(sum( matrix_3x3_of_neighbors)); % Count the matrix_3x3_of_neighbors + the central bin itself
        if ( sum_including_the_central_bin  > 0 ), % If the animal visited any of this bin's 8 neighbors (3x3 region)
            idx_timespent_density(ii_y_bin,ii_x_bin) = 1; % Put 1 in the central bin
        else  % If the animal did NOT visit any of this bin's 8 neighbors (3x3 region)
            idx_timespent_density(ii_y_bin,ii_x_bin) = 0; % Put 0 in the central bin (later we will divide by this 0 and will get NaN for the firing-rate map)
        end
    end
end
%                     warning on all

% Place Field = Spike Density / Time-Spent Density :
%                     warning off MATLAB:divideByZero ;
place_field_density_smoothed_with_NaN_X_Y = (place_field_density_smoothed_X_Y.* idx_timespent_density)./idx_timespent_density;
mat_timespent_density_smoothed_with_NaN_X_Y = (mat_timespent_density_smoothed_X_Y.* idx_timespent_density)./idx_timespent_density;
place_field_density_smoothed_with_NaN_normalized_X_Y = place_field_density_smoothed_with_NaN_X_Y./max(max(place_field_density_smoothed_with_NaN_X_Y));
[r,c] = size(place_field_density_smoothed_with_NaN_X_Y);
place_field_density_smoothed_with_NaN_5_level_binned_X_Y = ones(r,c)* NaN;
for ii = 1:r
    for jj = 1:c
        if ~isnan(place_field_density_smoothed_with_NaN_X_Y(ii,jj))
            position = histc(place_field_density_smoothed_with_NaN_normalized_X_Y(ii,jj),[0:0.2:1]);
            place_field_density_smoothed_with_NaN_5_level_binned_X_Y(ii,jj) = sum(position(1:end).*[0.0:0.2:1]);
        else
        end
    end
end

%                     warning on MATLAB:divideByZero ;

idx_notNaN_PlaceField_X_Y = find( ~isnan( place_field_density_smoothed_with_NaN_X_Y  ) ); % Find the indexes of non-NaN bins
idx_isNaN_PlaceField_X_Y = find( isnan( place_field_density_smoothed_with_NaN_X_Y  ) ); % Find the indexes of NaN bins

idx_notNaN_PlaceField_un_smoothed_rate_map_X_Y = find( ~isnan( place_field_density_raw_X_Y ) ); % Find the indexes of non-NaN bins

peak_firing_rate_X_Y = max(max(place_field_density_smoothed_with_NaN_X_Y));


% -----  INFORMATION PER SPIKE  (computed for the SMOOTHED field): -----
% Information_per_spike = sum( p_i * ( r_i / r ) * log2( r_i / r ) )
%    Where:
%       r_i = firing rate in bin i ;
%       p_i = occupancy of bin i = time-spent by bat in bin i / total time spent in all bins ;
%       r = mean( r_i ) = overall mean firing rate (mean over all the pixels)
% See: Skaggs WE, McNaughton BL, Wilson MA, Barnes CA, Hippocampus 6(2), 149-172 (1996).

r_i = place_field_density_smoothed_with_NaN_X_Y( idx_notNaN_PlaceField_X_Y ); % Use the SMOOTHED Place Field
p_i = mat_timespent_density_smoothed_with_NaN_X_Y( idx_notNaN_PlaceField_X_Y ) ./ ...
    sum( mat_timespent_density_smoothed_with_NaN_X_Y( idx_notNaN_PlaceField_X_Y ) ) ;
r_i = r_i(:) ; p_i = p_i(:) ; % Turn r_i and p_i into Column vectors
% % %             r = mean( r_i ) ;
r = sum( r_i .* p_i );
information_per_spike_X_Y = sum( p_i .* ( r_i / r ) .* log2( ( r_i + eps ) / r ) ) ; % I added a tiny number to avoid log(0)
SI_XY = information_per_spike_X_Y;
if figure_flag == 1
    figure;
    place_field_density_smoothed_with_NaN_X_Y_2 = imrotate(place_field_density_smoothed_with_NaN_X_Y,180);
    place_field_density_smoothed_with_NaN_X_Y_3 = fliplr(place_field_density_smoothed_with_NaN_X_Y_2);
    
    aa = imagesc(place_field_density_smoothed_with_NaN_X_Y_3);
    set(aa,'AlphaData',~isnan(place_field_density_smoothed_with_NaN_X_Y_3))
    axis tight; axis equal; axis off
    colormap(jet)
end

%% Rate map YZ
mat_spike_density_raw_Z_Y = zeros( diff(Z_room_dimensions)/bin_size_pixels_2_D,...
    diff(Y_room_dimensions)/bin_size_pixels_2_D ) + NaN ; % Initialize
mat_timespent_density_raw_Z_Y =  zeros( diff(Z_room_dimensions)/bin_size_pixels_2_D,...
    diff(Y_room_dimensions)/bin_size_pixels_2_D ) + NaN ; % Initialize
place_field_density_raw_Z_Y =  zeros( diff(Z_room_dimensions)/bin_size_pixels_2_D,...
    diff(Y_room_dimensions)/bin_size_pixels_2_D ) + NaN ; % Initialize

for ii_Y_bin = 1 : diff(Y_room_dimensions)/bin_size_pixels_2_D , % Loop over Y-bins
    for ii_z_bin = 1 : diff(Z_room_dimensions)/bin_size_pixels_2_D , % Loop over z-bins
        % Spike Density:
        mat_spike_density_raw_Z_Y( ii_z_bin,ii_Y_bin) = ... % All the data
            sum( y_spikes >= 1 + bin_size_pixels_2_D*(ii_Y_bin-1) & ...
            y_spikes <  1 + bin_size_pixels_2_D*(ii_Y_bin) & ...
            z_spikes >= 1 + bin_size_pixels_2_D*(ii_z_bin-1) & ...
            z_spikes <  1 + bin_size_pixels_2_D*(ii_z_bin)) ;
        % Time-Spent Density:
        mat_timespent_density_raw_Z_Y( ii_z_bin, ii_Y_bin ) = ...
            sum( y_video >= 1 + bin_size_pixels_2_D*(ii_Y_bin-1) & ...
            y_video <  1 + bin_size_pixels_2_D*(ii_Y_bin) & ...
            z_video >= 1 + bin_size_pixels_2_D*(ii_z_bin-1) & ...
            z_video <  1 + bin_size_pixels_2_D*(ii_z_bin) ) ;
        % Normalize Time-Spent Density from Video-Frames-Spent to Seconds-Spent
        mat_timespent_density_raw_Z_Y( ii_z_bin, ii_Y_bin ) = ...
            mat_timespent_density_raw_Z_Y( ii_z_bin, ii_Y_bin ) / frames_per_second ;
        %             Discard pixels in which the animal Spent less than a certain Minimal amount of time --
        %             (this is computed for the "idx_include_VT" data only, usually resulting in
        % DIFFERENT pixels being discarded for the Full data):
        if ( mat_timespent_density_raw_Z_Y( ii_z_bin, ii_Y_bin) < time_spent_minimum ),
            mat_timespent_density_raw_Z_Y( ii_z_bin, ii_Y_bin) = 0 ; % Discard this time-spent-density pixel
            mat_timespent_density_raw_Z_Y( ii_z_bin, ii_Y_bin ) = 0 ; % Discard this spike-density pixel
            
        end
        if(mat_timespent_density_raw_Z_Y~=0 & mat_timespent_density_raw_Z_Y ==0)
            disp(2);
        end
    end
end

% Place Field = Spike Density / Time-Spent Density :
%     warning off MATLAB:divideByZero ;
place_field_density_raw_Z_Y = mat_spike_density_raw_Z_Y ./ mat_timespent_density_raw_Z_Y;
%     warning on MATLAB:divideByZero ;


% Smoothing = convolve with gaussian kernel:
mat_spike_density_smoothed_Z_Y = imfilter(mat_spike_density_raw_Z_Y,gaussian_kernel);
mat_timespent_density_smoothed_Z_Y = imfilter(mat_timespent_density_raw_Z_Y,gaussian_kernel);


% Place Field smoothed = Spike Density smoothed / Time-Spent Density smoothed :
%     warning off MATLAB:divideByZero ;
place_field_density_smoothed_Z_Y = mat_spike_density_smoothed_Z_Y ./ mat_timespent_density_smoothed_Z_Y ;
%     warning on MATLAB:divideByZero ;

% ======= Compute the PF density with NaN's at unvisited location (will later be presented as white bins in the PF figure) : ==========

% "Legalize" a bin (remove NaN) if the bat visited any of the bin's 8 closest neighbours:
warning off all
idx_timespent_density = zeros( diff(Z_room_dimensions)/bin_size_pixels_2_D, diff(Y_room_dimensions)/bin_size_pixels_2_D ) + NaN ; % Initialize
for ii_Y_bin = 2 : (diff(Y_room_dimensions)/bin_size_pixels_2_D - 1) , % Loop over Y-bins, NOT INCL. CAMERA-VIEW EDGES
    for ii_z_bin = 2 : (diff(Z_room_dimensions)/bin_size_pixels_2_D - 1) , % Loop over z-bins, NOT INCL. CAMERA-VIEW EDGES
        matrix_3x3_of_neighbors = ...
            mat_timespent_density_raw_Z_Y( ii_z_bin-1 : ii_z_bin+1, ii_Y_bin-1 : ii_Y_bin+1 ) ;
        sum_including_the_central_bin = sum(sum( matrix_3x3_of_neighbors)); % Count the matrix_3x3_of_neighbors + the central bin itself
        if ( sum_including_the_central_bin  > 0 ), % If the animal visited any of this bin's 8 neighbors (3x3 region)
            idx_timespent_density(ii_z_bin,ii_Y_bin) = 1; % Put 1 in the central bin
        else  % If the animal did NOT visit any of this bin's 8 neighbors (3x3 region)
            idx_timespent_density(ii_z_bin,ii_Y_bin) = 0; % Put 0 in the central bin (later we will divide by this 0 and will get NaN for the firing-rate map)
        end
    end
end
warning on all

% Place Field = Spike Density / Time-Spent Density :
%     warning off MATLAB:divideByZero ;
place_field_density_smoothed_with_NaN_Z_Y = (place_field_density_smoothed_Z_Y.* idx_timespent_density)./idx_timespent_density;
mat_timespent_density_smoothed_with_NaN_Z_Y = (mat_timespent_density_smoothed_Z_Y.* idx_timespent_density)./idx_timespent_density;
place_field_density_smoothed_with_NaN_normalized_Z_Y = place_field_density_smoothed_with_NaN_Z_Y./max(max(place_field_density_smoothed_with_NaN_Z_Y));
[r,c] = size(place_field_density_smoothed_with_NaN_Z_Y);
place_field_density_smoothed_with_NaN_5_level_binned_Z_Y = ones(r,c)* NaN;
for ii = 1:r
    for jj = 1:c
        if ~isnan(place_field_density_smoothed_with_NaN_Z_Y(ii,jj))
            position = histc(place_field_density_smoothed_with_NaN_normalized_Z_Y(ii,jj),[0:0.2:1]);
            place_field_density_smoothed_with_NaN_5_level_binned_Z_Y(ii,jj) = sum(position(1:end).*[0.0:0.2:1]);
        else end
    end
end

%     warning on MATLAB:divideByZero ;

idx_notNaN_PlaceField_Z_Y = find( ~isnan( place_field_density_smoothed_with_NaN_Z_Y  ) ); % Find the indexes of non-NaN bins
idx_isNaN_PlaceField_Z_Y = find( isnan( place_field_density_smoothed_with_NaN_Z_Y  ) ); % Find the indexes of NaN bins

idx_notNaN_PlaceField_un_smoothed_rate_map_Z_Y = find( ~isnan( place_field_density_raw_Z_Y ) ); % Find the indexes of non-NaN bins

peak_firing_rate_Z_Y = max(max(place_field_density_smoothed_with_NaN_Z_Y));

% -----  INFORMATION PER SPIKE  (computed for the SMOOTHED field): -----
% Information_per_spike = sum( p_i * ( r_i / r ) * log2( r_i / r ) )
%    Where:
%       r_i = firing rate in bin i ;
%       p_i = occupancy of bin i = time-spent by bat in bin i / total time spent in all bins ;
%       r = mean( r_i ) = overall mean firing rate (mean over all the pixels)
% See: Skaggs WE, McNaughton BL, Wilson MA, Barnes CA, Hippocampus 6(2), 149-172 (1996).

r_i = place_field_density_smoothed_with_NaN_Z_Y( idx_notNaN_PlaceField_Z_Y ); % Use the SMOOTHED Place Field
p_i = mat_timespent_density_smoothed_with_NaN_Z_Y( idx_notNaN_PlaceField_Z_Y ) ./ ...
    sum( mat_timespent_density_smoothed_with_NaN_Z_Y( idx_notNaN_PlaceField_Z_Y ) ) ;
r_i = r_i(:) ; p_i = p_i(:) ; % Turn r_i and p_i into Column vectors
% % %             r = mean( r_i ) ;
r = sum( r_i .* p_i );
information_per_spike_Z_Y = sum( p_i .* ( r_i / r ) .* log2( ( r_i + eps ) / r ) ) ; % I added a tiny number to avoid log(0)

SI_YZ = information_per_spike_Z_Y;
if figure_flag == 1
    figure;
    place_field_density_smoothed_with_NaN_Z_Y_2 = imrotate(place_field_density_smoothed_with_NaN_Z_Y,180);
    place_field_density_smoothed_with_NaN_Z_Y_3 = fliplr(place_field_density_smoothed_with_NaN_Z_Y_2);
    aa = imagesc(place_field_density_smoothed_with_NaN_Z_Y_3);
    set(aa,'AlphaData',~isnan(place_field_density_smoothed_with_NaN_Z_Y_3))
    axis tight; axis equal; axis off
    colormap(jet)
end

%% XZ rate map
mat_spike_density_raw_X_Z = zeros( diff(Z_room_dimensions)/bin_size_pixels_2_D,...
    diff(X_room_dimensions)/bin_size_pixels_2_D ) + NaN ; % Initialize
mat_timespent_density_raw_X_Z =  zeros( diff(Z_room_dimensions)/bin_size_pixels_2_D,...
    diff(X_room_dimensions)/bin_size_pixels_2_D ) + NaN ; % Initialize
place_field_density_raw_X_Z =  zeros( diff(Z_room_dimensions)/bin_size_pixels_2_D,...
    diff(X_room_dimensions)/bin_size_pixels_2_D ) + NaN ; % Initialize

for ii_x_bin = 1 : diff(X_room_dimensions)/bin_size_pixels_2_D , % Loop over x-bins
    for ii_z_bin = 1 : diff(Z_room_dimensions)/bin_size_pixels_2_D , % Loop over z-bins
        % Spike Density:
        mat_spike_density_raw_X_Z( ii_z_bin,ii_x_bin) = ... % All the data
            sum( x_spikes >= 1 + bin_size_pixels_2_D*(ii_x_bin-1) & ...
            x_spikes <  1 + bin_size_pixels_2_D*(ii_x_bin) & ...
            z_spikes >= 1 + bin_size_pixels_2_D*(ii_z_bin-1) & ...
            z_spikes <  1 + bin_size_pixels_2_D*(ii_z_bin)) ;
        % Time-Spent Density:
        mat_timespent_density_raw_X_Z( ii_z_bin, ii_x_bin ) = ...
            sum( x_video >= 1 + bin_size_pixels_2_D*(ii_x_bin-1) & ...
            x_video <  1 + bin_size_pixels_2_D*(ii_x_bin) & ...
            z_video >= 1 + bin_size_pixels_2_D*(ii_z_bin-1) & ...
            z_video <  1 + bin_size_pixels_2_D*(ii_z_bin) ) ;
        % Normalize Time-Spent Density from Video-Frames-Spent to Seconds-Spent
        mat_timespent_density_raw_X_Z( ii_z_bin, ii_x_bin ) = ...
            mat_timespent_density_raw_X_Z( ii_z_bin, ii_x_bin ) / frames_per_second ;
        %             Discard pixels in which the animal Spent less than a certain Minimal amount of time --
        %             (this is computed for the "idx_include_VT" data only, usually resulting in
        % DIFFERENT pixels being discarded for the Full data):
        if ( mat_timespent_density_raw_X_Z( ii_z_bin, ii_x_bin) < time_spent_minimum ),
            mat_timespent_density_raw_X_Z( ii_z_bin, ii_x_bin) = 0 ; % Discard this time-spent-density pixel
            mat_timespent_density_raw_X_Z( ii_z_bin, ii_x_bin ) = 0 ; % Discard this spike-density pixel
            
        end
        if(mat_timespent_density_raw_X_Z~=0 & mat_timespent_density_raw_X_Z ==0)
            disp(2);
        end
    end
end

% Place Field = Spike Density / Time-Spent Density :
warning off MATLAB:divideByZero ;
place_field_density_raw_X_Z = mat_spike_density_raw_X_Z ./ mat_timespent_density_raw_X_Z;
warning on MATLAB:divideByZero ;


% Smoothing = convolve with gaussian kernel:
mat_spike_density_smoothed_X_Z = imfilter(mat_spike_density_raw_X_Z,gaussian_kernel);
mat_timespent_density_smoothed_X_Z = imfilter(mat_timespent_density_raw_X_Z,gaussian_kernel);


% Place Field smoothed = Spike Density smoothed / Time-Spent Density smoothed :
warning off MATLAB:divideByZero ;
place_field_density_smoothed_X_Z = mat_spike_density_smoothed_X_Z ./ mat_timespent_density_smoothed_X_Z ;
warning on MATLAB:divideByZero ;

% ======= Compute the PF density with NaN's at unvisited location (will later be presented as white bins in the PF figure) : ==========

% "Legalize" a bin (remove NaN) if the bat visited any of the bin's 8 closest neighbours:
warning off all
idx_timespent_density = zeros( diff(Z_room_dimensions)/bin_size_pixels_2_D, diff(X_room_dimensions)/bin_size_pixels_2_D ) + NaN ; % Initialize
for ii_x_bin = 2 : (diff(X_room_dimensions)/bin_size_pixels_2_D - 1) , % Loop over x-bins, NOT INCL. CAMERA-VIEW EDGES
    for ii_z_bin = 2 : (diff(Z_room_dimensions)/bin_size_pixels_2_D - 1) , % Loop over z-bins, NOT INCL. CAMERA-VIEW EDGES
        matrix_3x3_of_neighbors = ...
            mat_timespent_density_raw_X_Z( ii_z_bin-1 : ii_z_bin+1, ii_x_bin-1 : ii_x_bin+1 ) ;
        sum_including_the_central_bin = sum(sum( matrix_3x3_of_neighbors)); % Count the matrix_3x3_of_neighbors + the central bin itself
        if ( sum_including_the_central_bin  > 0 ), % If the animal visited any of this bin's 8 neighbors (3x3 region)
            idx_timespent_density(ii_z_bin,ii_x_bin) = 1; % Put 1 in the central bin
        else  % If the animal did NOT visit any of this bin's 8 neighbors (3x3 region)
            idx_timespent_density(ii_z_bin,ii_x_bin) = 0; % Put 0 in the central bin (later we will divide by this 0 and will get NaN for the firing-rate map)
        end
    end
end
warning on all

% Place Field = Spike Density / Time-Spent Density :
warning off MATLAB:divideByZero ;
place_field_density_smoothed_with_NaN_X_Z = (place_field_density_smoothed_X_Z.* idx_timespent_density)./idx_timespent_density;
mat_timespent_density_smoothed_with_NaN_X_Z = (mat_timespent_density_smoothed_X_Z.* idx_timespent_density)./idx_timespent_density;
place_field_density_smoothed_with_NaN_normalized_X_Z = place_field_density_smoothed_with_NaN_X_Z./max(max(place_field_density_smoothed_with_NaN_X_Z));
[r,c] = size(place_field_density_smoothed_with_NaN_X_Z);
place_field_density_smoothed_with_NaN_5_level_binned_X_Z = ones(r,c)* NaN;
for ii = 1:r
    for jj = 1:c
        if ~isnan(place_field_density_smoothed_with_NaN_X_Z(ii,jj))
            position = histc(place_field_density_smoothed_with_NaN_normalized_X_Z(ii,jj),[0:0.2:1]);
            place_field_density_smoothed_with_NaN_5_level_binned_X_Z(ii,jj) = sum(position(1:end).*[0.0:0.2:1]);
        else end
    end
end

warning on MATLAB:divideByZero ;

idx_notNaN_PlaceField_X_Z = find( ~isnan( place_field_density_smoothed_with_NaN_X_Z  ) ); % Find the indexes of non-NaN bins
idx_isNaN_PlaceField_X_Z = find( isnan( place_field_density_smoothed_with_NaN_X_Z  ) ); % Find the indexes of NaN bins

idx_notNaN_PlaceField_un_smoothed_rate_map_X_Z = find( ~isnan( place_field_density_raw_X_Z ) ); % Find the indexes of non-NaN bins

peak_firing_rate_X_Z = max(max(place_field_density_smoothed_with_NaN_X_Z));

% -----  INFORMATION PER SPIKE  (computed for the SMOOTHED field): -----
% Information_per_spike = sum( p_i * ( r_i / r ) * log2( r_i / r ) )
%    Where:
%       r_i = firing rate in bin i ;
%       p_i = occupancy of bin i = time-spent by bat in bin i / total time spent in all bins ;
%       r = mean( r_i ) = overall mean firing rate (mean over all the pixels)
% See: Skaggs WE, McNaughton BL, Wilson MA, Barnes CA, Hippocampus 6(2), 149-172 (1996).

r_i = place_field_density_smoothed_with_NaN_X_Z( idx_notNaN_PlaceField_X_Z ); % Use the SMOOTHED Place Field
p_i = mat_timespent_density_smoothed_with_NaN_X_Z( idx_notNaN_PlaceField_X_Z ) ./ ...
    sum( mat_timespent_density_smoothed_with_NaN_X_Z( idx_notNaN_PlaceField_X_Z ) ) ;
r_i = r_i(:) ; p_i = p_i(:) ; % Turn r_i and p_i into Column vectors
% % %             r = mean( r_i ) ;
r = sum( r_i .* p_i );
information_per_spike_X_Z = sum( p_i .* ( r_i / r ) .* log2( ( r_i + eps ) / r ) ) ; % I added a tiny number to avoid log(0)
SI_XZ = information_per_spike_X_Z;
if figure_flag == 1
    figure;
    place_field_density_smoothed_with_NaN_X_Z_2 = imrotate(place_field_density_smoothed_with_NaN_X_Z,180);
    place_field_density_smoothed_with_NaN_X_Z_3 = fliplr(place_field_density_smoothed_with_NaN_X_Z_2);
    aa = imagesc(place_field_density_smoothed_with_NaN_X_Z_3);
    set(aa,'AlphaData',~isnan(place_field_density_smoothed_with_NaN_X_Z_3))
    axis tight; axis equal; axis off
    colormap(jet)
end