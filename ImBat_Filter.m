function  [out_mov] = ImBat_Filter(mov);
% 
b = 10; %(smoothing )
% FiltB = mat2gray(FiltA-FBM);
% Scale data:


disp('scaling data');
mov = mat2gray(mov)*256;

disp(['smoothing data by a factor of ', num2str(b)]);
out_mov = (convn(mov, single(reshape([1 1 1] / b, 1, 1, [])), 'same'));

