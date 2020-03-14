function  [out_mov] = ImBat_Filter(mov,smth_frames);
%
b = 1; %(smoothing )
% FiltB = mat2gray(FiltA-FBM);
% Scale data:


disp('scaling data');

disp(['smoothing data by a factor of ', num2str(b)]);
% out_mov = medfilt3(mov,[1 1 b]);
out_mov = medfilt3(single(mov) ,[b b smth_frames]);

% spatial filter;


%out_mov = (convn(mov, single(reshape([1 1 1] / b, 1, 1, [])), 'same'));
