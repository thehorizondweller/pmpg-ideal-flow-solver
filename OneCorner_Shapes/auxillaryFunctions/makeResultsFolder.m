function makeResultsFolder()

    folderName = 'Results';
    if ~isfolder(folderName)
        mkdir(folderName);
        disp(['Folder "', folderName, '" created.']);
    else
        disp(['Folder "', folderName, '" exists. Clearing its contents...']);
        files = dir(fullfile(folderName, '*')); 
        for k = 1:length(files)
            if ~files(k).isdir 
                delete(fullfile(folderName, files(k).name)); 
            else
                if ~strcmp(files(k).name, '.') && ~strcmp(files(k).name, '..')
                    rmdir(fullfile(folderName, files(k).name), 's');
                end
            end
        end
    end
end