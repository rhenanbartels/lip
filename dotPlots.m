function [fig, xValues] = dotPlots(fig, data, offset, plotStyle,...
    markerFaceColor)

    if nargin == 3
        plotStyle = 'ko';
        markerFaceColor = 'k';
    elseif nargin == 4
        markerFaceColor = 'k';
            
    end   
          
    xValues = jitter(data, offset);
    figure(fig)
    
    hold on
    
    plot(xValues, data, plotStyle, 'MarkerFaceColor', markerFaceColor)
    
    hold off

end

function jitterAxisValues = jitter(data, offset)
    nData = length(data);
    jitterAxisValues = offset + 2 * randn(1, nData);
end