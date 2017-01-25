function [massVector, volumeVector] = volumeMassCurves(matFileNames,...
    lungPortion)

    if nargin == 1
        lungPortion = 'wholeLung';
    end
    
    if ~(strcmp(lungPortion, 'wholeLung') ||...
            strcmp(lungPortion, 'baseLung') ||...
            strcmp(lungPortion, 'middleLung') ||...
            strcmp(lungPortion, 'topLung'))
        error('Wrong lung portion');
    end

    nFiles = length(matFileNames);

    huValues = -1000:50;
    nHuValues = length(huValues);
    
    massVector = zeros(nFiles, nHuValues) * NaN;
    volumeVector = zeros(nFiles, nHuValues) * NaN;

    for i = 1:nFiles;
        currentMatFile = load(matFileNames{i});
        currentMatFile = getfield(currentMatFile, 'allResults');
        currentMatFile = getfield(currentMatFile, lungPortion);
        massVector(i, 1:size(currentMatFile.percentual.massPerDensity, 2))...
            = cumsum(currentMatFile.percentual.massPerDensity);
        volumeVector(i, 1:size(currentMatFile.percentual.volumePerDensity, 2))...
            = cumsum(...
            currentMatFile.percentual.volumePerDensity) /...
            sum(currentMatFile.percentual.volumePerDensity) * 100;
    end
end