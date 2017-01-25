function fileNames =  getMatFilesFullPath()
    pathName = uigetdir('Select the folder with .mat files');
    
    
    counter = 1;
    if pathName
        matFiles = dir([pathName filesep '*.mat']);
        nFiles = length(matFiles);
        for i = 1:nFiles
            currentFileName = matFiles(i).name;
            if ~startsWith(currentFileName, '._')
               fileNames{counter} = [pathName filesep currentFileName]; 
               counter = counter + 1;
            end
        end
    end
end


function match = startsWith(str, prefix)
    match = 0;
    
    n = length(prefix);
    if strcmp(str(1:n), prefix)
        match = 1;
    end
end