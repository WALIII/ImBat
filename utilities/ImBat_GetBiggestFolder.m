function out = ImBat_GetBiggestFolder(string);


OGDir = pwd;

    % Get all folders in directory
files2 = dir(pwd);
files2(ismember( {files2.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed

% Get a logical vector that tells which is a directory.
dirFlags2 = [files2.isdir];
% Extract only those that are directories.
subFolders2 = files2(dirFlags2);
% Print folder names to command window.
for k = 1 : length(subFolders2)
	fprintf('Sub folder #%d = %s\n', k, subFolders2(k).name);
    folds{k} = subFolders2(k).name;
end

u = regexp(folds,string);


    if size(find([u{:}]>0),2) >0;
         flys = find([u{:}]>0); % these are the 'flight folders'
    else
        flys = 1:length(subFolders2); % Take all flights
    end
    
    for ii = 1:size(flys,2);
        cd(subFolders2(flys(ii)).name);
        indX(ii) = DirSize(pwd);
        cd(OGDir);
    end
    
    [a b] = sort(indX,'descend');
    
    out = subFolders2(flys(b)).name;
    
    % get size of each directry 
    
    
    % get the biggest one:
    
    