function fig = plotStatsLines(fig, averageMass, lowerStd, upperStd, statsColor)
    
    if nargin == 4
        statsColor = [0.5, 0.5, 0.5];    
    end
    
    figure(fig)
    hold on
    
    lineWidth = 1.2;
    plot(averageMass, 1:100, 'Color', statsColor, 'LineWidth', lineWidth)
    plot(lowerStd, 1:100, 'Color', statsColor, 'LineWidth', lineWidth)
    plot(upperStd, 1:100, 'Color', statsColor, 'LineWidth', lineWidth)
    
    hold off
end
