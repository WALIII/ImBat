function [out_mov] = ImBat_TemporalDownSample(mov,ds_temp);
% grouped z project


LastFrame = size(mov,3);

frameIdx = 1:ds_temp:LastFrame;
counter = 1;
mov = (convn(mov, single(reshape([1 1 1] / ds_temp, 1, 1, [])), 'same'));

for i = 1:size(frameIdx,2)-1
  out_mov(:,:,counter) = mean(mov(:,:,frameIdx(i):frameIdx(i+1)),3);
counter= counter+1;
end
