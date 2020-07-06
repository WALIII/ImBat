function A = SimpleRead

fname = 'G.tiff';
info = imfinfo(fname);
num_images = numel(info);
for k = 1:num_images
    A(:,:,k) = single(imread(fname, k, 'Info', info));
end