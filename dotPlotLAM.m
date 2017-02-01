function dotPlotLAM()
    matFileNames = getMatFilesFullPath();
    [massVector, volumeVector] = volumeMassCurves(matFileNames);
    
    [m3, m15, m85, m97] = calculateMs(massVector, volumeVector);
    
    fig = figure;
    
    [fig, xValuesM15] = dotPlots(fig, m15, 10, 'ko', 'k');
    [fig, xValuesM3] = dotPlots(fig, m3, 20, 'ko', 'w');
    
    fig = averageLineDotPlot(fig, m15, xValuesM15);
    fig = averageLineDotPlot(fig, m3, xValuesM3);
    
    set(gca, 'xtick', [10, 20]);
    set(gca, 'xticklabels', {'LAM', 'Control'})
end