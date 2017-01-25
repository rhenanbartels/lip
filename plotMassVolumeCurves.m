function fig = plotMassVolumeCurves(fig, massVector, volumeVector,...
    plotColor)

    if nargin == 3
        plotColor = 'k';
    end

    figure(fig);
    hold on
    
    for i = 1:size(massVector, 1)
     plot(massVector(i, ~isnan(massVector(i, :))),...
         volumeVector(i, ~isnan(volumeVector(i, :))), plotColor)
    end
    
    hold off

end