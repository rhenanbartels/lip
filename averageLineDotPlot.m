function fig = averageLineDotPlot(fig, data, xValues, plotColor)

    if nargin == 3
        plotColor = 'r';
    end
    
    figure(fig)
    
    hold on
    plot([min(xValues), max(xValues)], [mean(data), mean(data)], plotColor)
    hold off
end