function [massVector, volumeVector] = volumeMassCurves(matFileNames)

    nFiles = length(matFileNames);

    huValues = -1000:50;
    nHuValues = length(huValues);
    
    massVector = zeros(nFiles, nHuValues);
    volumeVector = zeros(nFiles, nHuValues);

    for i = 1:nFiles;
        currentMatFile = load(matFileNames{i});
        currentMatFile = getfield(currentMatFile, 'allResults');
        massVector(i, :) = cumsum(currentMatFile.wholeLung.percentual.massPerDensity);
        volumeVector(i, :) = cumsum(...
            currentMatFile.wholeLung.percentual.volumePerDensity) /...
            sum(currentMatFile.wholeLung.percentual.volumePerDensity) * 100;
    end
end