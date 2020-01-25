function  [out_mov] = ImBat_Filter(mov,smth_frames);
%
b = 3; %(smoothing )
% FiltB = mat2gray(FiltA-FBM);
% Scale data:


disp('scaling data');

disp(['smoothing data by a factor of ', num2str(b)]);
% out_mov = medfilt3(mov,[1 1 b]);
out_mov = medfilt3(single(mov) ,[1 1 smth_frames]);

out_mov  = round(mat2gray(out_mov))*256;

%out_mov = (convn(mov, single(reshape([1 1 1] / b, 1, 1, [])), 'same'));
