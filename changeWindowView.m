function [displayLow, displayHigh] = changeWindowView(dicomData, wl, ww)

    if nargin < 2
        ww = 1400;
        wl = -500;
    elseif nargin < 3
        ww = 1400;        
    end
    
    displayLow = max(wl - 0.5 * ww, min(dicomData(:)));
    displayHigh = max(wl + 0.5 * ww, min(dicomData(:)));
    
    %set(axesObj, 'CLim', [displayLow, displayHigh]);
end