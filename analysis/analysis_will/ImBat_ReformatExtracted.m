function ImBat_ReformatExtracted
% ImBat_ReformatExtracted.m

% hopefully, one time use function to reformat the processed data folders

% WAL3
% d12/12/2020


% check if any mov has been extracted,if not,
if exist('extracted')>1
    
else % move files into extracted folder
    mkdir('extracted');
    A= [subFolders(i).folder,'/',subFolders(i).name]
    B= [A,'/extracted']
    f_raw=dir(fullfile(A));
    for i = 3:size(f_raw,1)-1
        movefile(fullfile(f_raw(i).folder,f_raw(i).name),B);
    end
end
