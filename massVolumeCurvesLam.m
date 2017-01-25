function massVolumeCurvesLam()
    matFileNames = getMatFilesFullPath();
    
    [massVector, volumeVector] = volumeMassCurves(matFileNames);
    [averageMass, lowerStd, upperStd] = calculateAverageAndStd(massVector,...
        volumeVector);
    
    fig = figure;
    
    fig = plotMassVolumeCurves(fig, massVector, volumeVector);
    fig = plotStatsLines(fig, averageMass, lowerStd, upperStd);
    
    xlabel('Mass (g)');
    ylabel('Volume (%)');
end