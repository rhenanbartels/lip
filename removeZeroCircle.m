function roi = removeZeroCircle(roi)    
    ncols = size(roi, 2);
    
    for j = 1:ncols
        roi(:, j) = discoverBoundary(roi(:, j));
    end
    
end

function roiColumn = discoverBoundary(roiColumn)
    nrows = length(roiColumn);
    upper = false;
    lower = false;
    
    upperBoundary = 1;
    lowerBoundary = nrows;
    
    i = 1;
    while ~(upper && lower)
        if roiColumn(i) ~= 0 && ~upper
           upper = true; 
           upperBoundary = i;
        end
        if roiColumn(end - i) ~= 0 && ~lower
            lower = true;
            lowerBoundary = nrows - i;
        end
        i = i + 1;
        if i == nrows
            upper = true;
            lower = true;
        end
    end
    if upperBoundary == 1 && lowerBoundary == nrows
        roiColumn = ones(size(roiColumn)) * NaN;
    else
        roiColumn(1:upperBoundary - 1) = NaN;
        roiColumn(lowerBoundary + 1:end) = NaN;
    end
end