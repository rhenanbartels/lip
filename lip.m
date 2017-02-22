function lip()
    drawInterface()
end

function [metadata, dicomImages] = getDicomData(dirName)
    filesAndFolders = dir(dirName);
    nItems = length(filesAndFolders);
    
    
   % h = waitbar(1,'Loading Metadata and DICOMS...');
    
    counter = 1;
    for idx = 1:nItems
        currentItemName = filesAndFolders(idx).name;
        if ~strcmp(currentItemName, '.') && ~strcmp(currentItemName, '..')...
                && ~strcmp(currentItemName, '.DS_Store')
            try
                fullPath = [dirName filesep currentItemName];
                metadata{counter} = dicominfo(fullPath);
                dicomImages(:, :, counter) = dicomread(fullPath);
                counter = counter + 1;
                %waitbar(counter / nItems);
            catch
                continue
            end
        end
    end
    %close(h)
end

%%% DIVIDE LUNG %%%
function [firstLevel, secondLevel, thirdLevel] = divideLungIndexes(mask)
    %numberOfMasks = unique(mask);    
    
    %if length(numberOfMasks) == 2
        nSlices = size(mask, 3);
        [firstLevel, secondLevel, thirdLevel] = lungThresholds(nSlices);
    %end
end

function [firstLevel, secondLevel, thirdLevel] = lungThresholds(nSlices)
    %Centimeters lung size;
    oneThird = nSlices / 3;
    floorOneThird = floor(oneThird);

    diffOneThird = abs(floorOneThird - oneThird);

    if (diffOneThird)
        firstLevel = [1, floorOneThird + 1];
        if diffOneThird > 0 && diffOneThird < 0.5
            secondLevel = [firstLevel(2) + 1, firstLevel(2) + floorOneThird];
        else
            secondLevel = [firstLevel(2) + 1, firstLevel(2) + floorOneThird + 1];
        end

    else
        firstLevel = [1 floorOneThird];
        secondLevel = [floorOneThird + 1, firstLevel(2) + floorOneThird];
    end

    thirdLevel = [secondLevel(2) + 1, secondLevel(2) + floorOneThird];
    assert(thirdLevel(2) == nSlices, 'Foi diferente')
end

%%% CALCULATIONS %%%
function [resultsClassicalWholeLung, resultsPercentileWholeLung,...
    resultsClassicalBaseLung, resultsPercentileBaseLung,...
    resultsClassicalMiddleLung, resultsPercentileMiddleLung,...
    resultsClassicalTopLung, resultsPercentileTopLung] =...
    allAnalysis(lung, lungMask, metadata)

   %Discover where the mask
    maskPosition = find(sum(sum(lungMask)) >= 1);
    mask = lungMask(:, :, maskPosition);
    lung = lung(:, :, maskPosition);
    
    %Divide Lung indexes
    [firstLevel, secondLevel, thirdLevel] = divideLungIndexes(mask);
    
    %Whole lung
    resultsClassicalWholeLung = classicalAnalysis(lung, mask, metadata);
    resultsPercentileWholeLung = lungAnalysis(lung, mask, metadata);
    
    %Base Lung
    resultsClassicalBaseLung = classicalAnalysis(lung(:, :,...
        firstLevel(1):firstLevel(2)), mask(:, :,...
        firstLevel(1):firstLevel(2)), metadata);
    resultsPercentileBaseLung = lungAnalysis(lung(:, :,...
        firstLevel(1):firstLevel(2)), mask(:, :,...
        firstLevel(1):firstLevel(2)), metadata);    
    
    %Middle Lung
    resultsClassicalMiddleLung = classicalAnalysis(lung(:, :,...
        secondLevel(1):secondLevel(2)), mask(:, :,...
        secondLevel(1):secondLevel(2)), metadata);
    resultsPercentileMiddleLung = lungAnalysis(lung(:, :,...
        secondLevel(1):secondLevel(2)), mask(:, :,...
        secondLevel(1):secondLevel(2)), metadata);
    
    %Top Lung
    resultsClassicalTopLung = classicalAnalysis(lung(:, :,...
        thirdLevel(1):thirdLevel(2)), mask(:, :,...
        thirdLevel(1):thirdLevel(2)), metadata);
    resultsPercentileTopLung = lungAnalysis(lung(:, :,...
        thirdLevel(1):thirdLevel(2)), mask(:, :,...
        thirdLevel(1):thirdLevel(2)), metadata);
end

function [lungOrMask, metadatas] = sortLungOrMask(lungOrMask, metadatas)
    nSlices = length(metadatas);
    sliceLocations = zeros(1, nSlices);
    
    for idx = 1:nSlices
        sliceLocations(idx) = metadatas{idx}.SliceLocation;
    end
    
    [trash, sortedIndexes] = sort(sliceLocations);
    
%    lungOrMask = lungOrMask(:, :, sortedIndexes);
    if nargout == 2
        metadatas = metadatas(sortedIndexes);
    end
end

function results  = lungAnalysis(lung, mask, metadata)

    voxelVolume = calculateVoxelVolume(metadata{1}, metadata{2});

    roiAir = -1000;
    roiTissue = 50;
    
    lung(mask == 0) = 10000;
    
    lung = int16(lung);

    lung(lung < -1000 | lung > 50) = [];



    huValues = double(unique(lung));
    nLung = size(lung, 3);

    massPerDensity = zeros(1, length(huValues));
    volumePerDensity = zeros(1, length(huValues));
    voxelPerDensity = zeros(1, length(huValues));
    counter = 1;    

    %h = waitbar(0, 'Calculating Mass Percentile...');

    for hu = huValues    
        nVoxels = length(lung(lung == hu));
        massPerDensity(counter) = ((hu - roiAir)/(roiTissue - roiAir))...
            * voxelVolume * 1.04 * nVoxels;
        volumePerDensity(counter) = nVoxels * voxelVolume / 1000;
        voxelPerDensity(counter) = nVoxels;
        lung(lung == hu) = [];
        counter = counter + 1;
        %waitbar(counter / length(huValues))
    end
   % close(h)
    sumMass = cumsum(massPerDensity);

    sumVolume = cumsum(volumePerDensity);
    cumulativeVolume = sumVolume / sumVolume(end) * 100;
    %cumulativeVoxel = cumsum(voxelPerDensity);

    [~, pos3] = min(abs(cumulativeVolume - 3));
    [~, pos15] = min(abs(cumulativeVolume - 15));    
    [~, pos65] = min(abs(cumulativeVolume - 65));    
    [~, pos70] = min(abs(cumulativeVolume - 70));    
    [~, pos75] = min(abs(cumulativeVolume - 75));    
    [~, pos80] = min(abs(cumulativeVolume - 80));
    [~, pos85] = min(abs(cumulativeVolume - 85));
    [~, pos97] = min(abs(cumulativeVolume - 97));

    results.p3Mass = sumMass(pos3);
    results.p15Mass = sumMass(pos15);
    results.p65Mass = sumMass(end) - sumMass(pos65);
    results.p70Mass = sumMass(end) - sumMass(pos70);
    results.p75Mass = sumMass(end) - sumMass(pos75);
    results.p80Mass = sumMass(end) - sumMass(pos80);
    results.p85Mass = sumMass(end) - sumMass(pos85);
    results.p97Mass = sumMass(end) - sumMass(pos97);

    results.p65Vol = sumVolume(end) - sumVolume(pos65);
    results.p70Vol = sumVolume(end) - sumVolume(pos70);
    results.p75Vol = sumVolume(end) - sumVolume(pos75);
    results.p80Vol = sumVolume(end) - sumVolume(pos80);

    results.p3Volume = sumVolume(pos3);
    results.p15Volume = sumVolume(pos15);
    results.p85Volume = sumVolume(end) - sumVolume(pos85);
    results.p97Volume = sumVolume(end) - sumVolume(pos97);
    
    results.p3Density = results.p3Mass / results.p3Volume;
    results.p15Density = results.p15Mass / results.p15Volume;
    results.p85Density = results.p85Mass / results.p85Volume;
    results.p97Density = results.p97Mass /results. p97Volume;

    results.totalVolume = sumVolume(end);
    results.totalMass = sumMass(end);
    
    results.massPerDensity = massPerDensity;
    results.volumePerDensity = volumePerDensity;
    results.huValues = huValues;
    results.voxelPerDensity = voxelPerDensity;
end

function results = classicalAnalysis(lung, mask, metadata)
    voxelVolume = calculateVoxelVolume(metadata{1}, metadata{2}); 
    
    
    roiAir = -1000;
    roiTissue = 50;
    
    lung(mask == 0) = 10000;
    
    nSlices = size(lung, 3);
   
    
    hyper = [-1000, -901];
    normally = [-900, -501];
    poor = [-500, -101];
    non = [-100, 100];
    
    %Mass
    hyperMassPerSlice = zeros(1, nSlices);
    poorMassPerSlice = zeros(1, nSlices);
    normallyMassPerSlice = zeros(1, nSlices);
    nonMassPerSlice = zeros(1, nSlices);
    totalMassPerSlice = zeros(1, nSlices);
    
    hyperMassPercentualPerSlice = zeros(1, nSlices);
    poorMassPercentualPerSlice = zeros(1, nSlices);
    normallyMassPercentualPerSlice = zeros(1, nSlices);
    nonMassPercentualPerSlice = zeros(1, nSlices);
    
    %Volume
    hyperVolumePerSlice = zeros(1, nSlices);
    poorVolumePerSlice = zeros(1, nSlices);
    normallyVolumePerSlice = zeros(1, nSlices);
    nonVolumePerSlice = zeros(1, nSlices);
    totalVolumePerSlice = zeros(1, nSlices);
    
    hyperVolumePercentualPerSlice = zeros(1, nSlices);
    poorVolumePercentualPerSlice = zeros(1, nSlices);
    normallyVolumePercentualPerSlice = zeros(1, nSlices);
    nonVolumePercentualPerSlice = zeros(1, nSlices);
    
    %Density
    hyperDensityPerSlice = zeros(1, nSlices);
    poorDensityPerSlice = zeros(1, nSlices);
    normallyDensityPerSlice = zeros(1, nSlices);
    nonDensityPerSlice = zeros(1, nSlices);
    totalDensityPerSlice = zeros(1, nSlices);    
    
    %h = waitbar(0, 'Calculating Classical Indexes...');
    
    for idx = 1:nSlices
        currentSlice = lung(:, :, idx);
        currentMask = mask(:, :, idx) >= 1;
        voxels = single(int32(currentSlice(currentMask)));
        
        
        % Group the voxels according to their HU value.
        voxels_non  = voxels(voxels >= non(1) & voxels <=...
            non(2));
        voxels_poor  = voxels(voxels >= poor(1) & voxels <=...
            poor(2));
        voxels_norm  = voxels(voxels >= normally(1) & voxels <=...
            normally(2));
        voxels_hyp  = voxels(voxels >= hyper(1) & voxels <=...
            hyper(2));
        voxels_total = voxels(voxels >= hyper(1) & voxels <=...
            non(2));
        
        %Mass
        hyperMassPerSlice(idx) = sum(((voxels_hyp  - roiAir)/(roiTissue -...
            roiAir)) * voxelVolume * 1.04);
        
        poorMassPerSlice(idx) = sum(((voxels_poor  - roiAir)/(roiTissue -...
            roiAir)) * voxelVolume * 1.04);
        
        normallyMassPerSlice(idx) = sum(((voxels_norm  - roiAir)/(roiTissue -...
            roiAir)) * voxelVolume * 1.04);
        
        nonMassPerSlice(idx) = sum(((voxels_non  - roiAir)/(roiTissue -...
            roiAir)) * voxelVolume * 1.04);
        
        totalMassPerSlice(idx) = sum(((voxels_total  - roiAir)/(roiTissue -...
            roiAir)) * voxelVolume * 1.04);
        
        %Mass Percentuais
        hyperMassPercentualPerSlice(idx) = hyperMassPerSlice(idx) /...
            totalMassPerSlice(idx) * 100;
        
        poorMassPercentualPerSlice(idx) = poorMassPerSlice(idx) / ...
            totalMassPerSlice(idx) * 100;
        
        normallyMassPercentualPerSlice(idx) = normallyMassPerSlice(idx) / ...
            totalMassPerSlice(idx) * 100;
        
        nonMassPercentualPerSlice(idx) = nonMassPerSlice(idx) / ...
            totalMassPerSlice(idx) * 100;
        
        
         
        %Volume
        hyperVolumePerSlice(idx) = length(voxels_hyp) * voxelVolume;
        
        poorVolumePerSlice(idx) = length(voxels_poor) * voxelVolume;
        
        normallyVolumePerSlice(idx) = length(voxels_norm) * voxelVolume;
        
        nonVolumePerSlice(idx) = length(voxels_non) * voxelVolume;
        
        totalVolumePerSlice(idx) = length(voxels_total) * voxelVolume;
        
        %Volume Percentuais
        hyperVolumePercentualPerSlice(idx) = hyperVolumePerSlice(idx) /...
            totalVolumePerSlice(idx) * 100;
        
        poorVolumePercentualPerSlice(idx) = poorVolumePerSlice(idx) / ...
            totalVolumePerSlice(idx) * 100;
        
        normallyVolumePercentualPerSlice(idx) = normallyVolumePerSlice(idx) / ...
            totalVolumePerSlice(idx) * 100;
        
        nonVolumePercentualPerSlice(idx) = nonVolumePerSlice(idx) / ...
            totalVolumePerSlice(idx) * 100;
        
        
        %Density
        hyperDensityPerSlice(idx) = hyperMassPerSlice(idx) /...
            hyperVolumePerSlice(idx);
        
        poorDensityPerSlice(idx) = poorMassPerSlice(idx) /...
            poorVolumePerSlice(idx);
        
        normallyDensityPerSlice(idx) = normallyMassPerSlice(idx) /...
            normallyVolumePerSlice(idx);
        
        nonDensityPerSlice(idx) = nonMassPerSlice(idx) /...
            nonVolumePerSlice(idx);
        
        totalDensityPerSlice(idx) = totalMassPerSlice(idx) /...
            totalVolumePerSlice(idx);
        
        %waitbar(idx / nSlices)
    end    
    %close(h)
    
    %Whole Lung
    %Mass
    results.hyperMass = sum(hyperMassPerSlice);
    results.poorMass = sum(poorMassPerSlice);
    results.nonMass = sum(nonMassPerSlice);
    results.normallyMass = sum(normallyMassPerSlice);
    results.totalMass = sum(totalMassPerSlice);
    
    %Volume
    results.hyperVolume = sum(hyperVolumePerSlice);
    results.poorVolume = sum(poorVolumePerSlice);
    results.nonVolume = sum(nonVolumePerSlice);
    results.normallyVolume = sum(normallyVolumePerSlice);
    results.totalVolume = sum(totalVolumePerSlice);
    
    %Density
    results.hyperDensity = results.hyperMass / results.hyperVolume;
    results.poorDensity = results.poorMass / results.poorVolume;
    results.nonDensity = results.nonMass / results.nonVolume;
    results.normallyDensity = results.normallyMass / results.normallyVolume;
    results.totalDensity = results.totalMass / results.totalVolume;
    
    
    %Assign results
    %Mass
    results.hyperMassPerSlice = hyperMassPerSlice;
    results.poorMassPerSlice = poorMassPerSlice;
    results.normallyMassPerSlice = normallyMassPerSlice;
    results.nonMassPerSlice = nonMassPerSlice;
    results.totalMassPerSlice = totalMassPerSlice;
    
    results.hyperMassPercentualPerSlice = hyperMassPercentualPerSlice;
    results.poorMassPercentualPerSlice = poorMassPercentualPerSlice;
    results.normallyMassPercentualPerSlice = normallyMassPercentualPerSlice;
    results.nonMassPercentualPerSlice = nonMassPercentualPerSlice;
    
    %Volume
    results.hyperVolumePerSlice = hyperVolumePerSlice;
    results.poorVolumePerSlice = poorVolumePerSlice;
    results.normallyVolumePerSlice = normallyVolumePerSlice;
    results.nonVolumePerSlice = nonVolumePerSlice;
    results.totalVolumePerSlice = totalVolumePerSlice;
    
    results.hyperVolumePercentualPerSlice = hyperVolumePercentualPerSlice;
    results.poorVolumePercentualPerSlice = poorVolumePercentualPerSlice;
    results.normallyVolumePercentualPerSlice = normallyVolumePercentualPerSlice;
    results.nonVolumePercentualPerSlice = nonVolumePercentualPerSlice;
    
    %Density
    results.hyperDensityPerSlice = hyperDensityPerSlice;
    results.poorDensityPerSlice = poorDensityPerSlice;
    results.normallyDensityPerSlice = normallyDensityPerSlice;
    results.nonDensityPerSlice = nonDensityPerSlice;
    results.totalDensityPerSlice = totalDensityPerSlice;
    
end

function patientName = getPatientName(metadata)
    patientName = 'patient_name';
    if isfield(metadata, 'PatientName')
        if isfield(metadata.PatientName,'FamilyName')
            patientName = metadata.PatientName.FamilyName;
        end
    end
end

%%% Voxel Volume %%%
function voxelVolume = calculateVoxelVolume(metadata, metadata2)
if isfield(metadata,'SpacingBetweenSlices');
    if isfield(metadata,'SliceThickness')
        if abs(metadata.SpacingBetweenSlices) < metadata.SliceThickness
            voxelVolume =(metadata.PixelSpacing(1) ^ 2 *...
                metadata.SliceThickness * 0.001) *...
                (abs(metadata.SpacingBetweenSlices) / metadata.SliceThickness);
        else
            voxelVolume = (metadata.PixelSpacing(1) ^ 2 *...
                metadata.SliceThickness * 0.001);
        end
    else
        voxelVolume = (metadata.PixelSpacing(1) ^ 2 *...
            metadata.SliceThickness * 0.001);
    end
else
    
    if isfield(metadata,'SliceThickness')==1;
        thick=abs(metadata.SliceThickness);
    elseif isfield(metadata,'SpacingBetweenSlices');
        thick=abs(metadata.SpacingBetweenSlices);
    else
        thick=abs(metadata.PixelDimensions(3));
    end
    
    if isfield(metadata, 'SliceLocation')
        SpacingBetweenSlices = abs(metadata2.SliceLocation -...
            metadata.SliceLocation);
    end

    SliceThickness = metadata.SliceThickness;
    voxelVolume = (metadata.PixelSpacing(1) ^ 2 * thick * 0.001) *...
        (SpacingBetweenSlices / SliceThickness);
end
end

%%%% LUNG CALIBRATION %%%%
function calibratedLung = lungCalibration(rawLung, mAir, mTissue)
    coef = polyfit([round(mAir), round(mTissue)], [-1000 50], 1);
    calibratedLung = int16(rawLung) * int16(coef(1)) + int16(coef(2));

end

function masks = getHDRMask(fileName)
    masks = analyze75read(fileName);
end

function showDicom(axes_handle, lung_slice)
    axes(axes_handle);
    imagesc(lung_slice);
    colormap(gray)
    axis equal;axis off 
end

function prepareAndShowLung(handles, currentSlice)
    axes_handle = handles.gui.mainAxes;    
    lungSlice = handles.data.lung(:, :, currentSlice);
    showDicom(axes_handle, lungSlice)
    isToShowMask = get(handles.gui.showMask, 'Value');
    if isToShowMask
        lungDim = size(lungSlice, 1);
        mask = handles.data.lungMask(:, :, currentSlice);
        showMask(lungDim, mask)
    end
end

function showMask(lungDim, mask)
    mask = mask >= 1;
 
    color1 = 0; color2 = 0.8; color3 = 0;
    colorMask = cat(3, color1 * ones(lungDim), color2 * ones(lungDim),...
         color3 * ones(lungDim));
    
    hold on
    h = imshow(colorMask);
    set(h, 'AlphaData', mask);    
end

function currentSlicePosition = getCurrentSlicePosition(handles)
    currentSlicePosition = str2double(get(handles.gui.currentSlicePosition,...
        'String'));
end

function setCurrentSlicePosition(handles, newSlice)
    set(handles.gui.currentSlicePosition,...
     'String', num2str(newSlice));
end

function makeWidgetsVisible(handles)
    set(handles.gui.buttonNext, 'Visible', 'On')
    set(handles.gui.buttonPrevious, 'Visible', 'On')
    set(handles.gui.currentSlicePosition, 'Visible', 'On')
    set(handles.gui.patName, 'Visible', 'On') 
    set(handles.gui.maskMenu, 'Enable', 'On')    
end

function lung = uncalibrateLung(lung, metadata)
    lung = (single(lung) - metadata.RescaleIntercept) / metadata.RescaleSlope;
end


function [hasROIS, avgAir, avgTissue, roiAir, roiTissue] =...
    checkForCalibrationROI(handles)   
    
    avgAir = NaN; avgTissue = NaN; roiAir = NaN; roiTissue = NaN;

    children = get(handles.gui.mainAxes, 'Child');
    nCircles = 0;
    for i = 1:length(children)
        fields = get(children(i));
        if isfield(fields, 'Type') && strcmp(fields.Type, 'rectangle')  
            nCircles = nCircles + 1;            
        end
    end
    hasROIS = nCircles == 2;
    
    if hasROIS
        [avgAir, roiAir] = averageCircle(handles, 'air');
        [avgTissue, roiTissue] = averageCircle(handles, 'tissue');
    end
end

function [m, imgMask] = averageCircle(handles, roi_type)
%meanDisk computes mean of values inside a circle
%   M = meanDisk(IMG, XC, YC, R) returns the mean of IMG(Y,X) for all X and
%   Y such that the Euclidean distance of (X,Y) from (XC,YC) is less than
%   R. IMG must be 2-D, R must be positive, and some elements of IMG must
%   lie within the circle.
% This section is for efficiency only - avoids wasting computation time on
% pixels outside the bounding square

switch roi_type
    case 'air'
        circle_object = findobj(gca, 'Type', 'rectangle','-and',...
            'EdgeColor', 'g');
        roi_handle = handles.roi_air_properties;        
    otherwise
        circle_object = findobj(gca, 'Type', 'rectangle','-and',...
            'EdgeColor', 'r');
        roi_handle = handles.roi_tissue_properties;
end
position = get(circle_object, 'Position');
%x do circulo
r = position(3) / 2;
xc = position(1) + r; %reposiciona o centro do circulo.
%y do circulo
yc = position(2) + r; %reposiciona o centro do circulo.
%raio

img = handles.data.lung;

slice = str2double(get(handles.gui.currentSlicePosition, 'string'));
[sy sx] = size(img(:,:, slice));
xmin = max(1, floor(xc-r));
xmax = min(sx, ceil(xc+r));
ymin = max(1, floor(yc-r));
ymax = min(sy, ceil(yc+r));
img = img(ymin:ymax, xmin:xmax, slice); % trim boundaries
%figure;imagesc(img);colormap(gray(156))
xc = xc - xmin + 1;
yc = yc - ymin + 1;
% Make a circle mask
[x y] = meshgrid(1:size(img,2), 1:size(img,1));
mask = (x-xc).^2 + (y-yc).^2 < r.^2;
% Compute mean
m = sum(sum(double(img) .* mask)) / sum(mask(:));
imgMask = double(img) .* mask;
end

function [a, properties] = drag_and_drop(x, y, r, color)
d = r*2;
px = x;
py = y;
dragging = [];
orPos = [];

set(gcf,'WindowButtonUpFcn',@dropObject,'units','normalized',...
    'WindowButtonMotionFcn',@moveObject);

a = rectangle('Position',[px py d d],'Curvature',[1,1],...
    'ButtonDownFcn',@dragObject,'EdgeColor', color, 'LineWidth', 2);


properties = get(a, 'Position');

    function dragObject(hObject,eventdata)
        selection_type = get(gcf, 'SelectionType');
        dragging = hObject;
        orPos = get(a,'Position');
        face_color = get(a, 'FaceColor');
        if strcmp(selection_type, 'open')
            if isempty(face_color) || strcmp('none', face_color)
                set(a, 'FaceColor', color)                
            else
                set(a, 'FaceColor', 'None')               
            end
        end
        
    end
    function dropObject(hObject,eventdata)
        
        dragging = [];
        handles = guidata(hObject);
        color = get(a, 'EdgeColor');
        if color == [0 1 0];
            [m, imgMask] = averageCircle(handles, 'air');
            disp(['Valor médio da ROI da TRAQUEIA é: ' num2str(round(m))]);
        else
            [m, imgMask] = averageCircle(handles, 'tissue');
            disp(['Valor médio da ROI da AORTA é: ' num2str(round(m))]);
        end
        
    end
    function moveObject(hObject,eventdata)
        selection_type = get(gcf,'SelectionType');
        newPos = get(gca,'CurrentPoint');
        %current_position = get(a, 'Position');        
        if ~isempty(dragging)
            if ~strcmp('extend', selection_type)
                current_position = get(a, 'Position');
                %Avoid to the circle to flip to the right.
                posDiff(1) = newPos(1, 1) - orPos(1, 1);
                posDiff(2) = -(newPos(1, 2) - orPos(1, 2));
                current_position = current_position +...
                    [[posDiff(1) -posDiff(2)] 0 0];
                if d > 0
                    %current_position(3:4) = d;
                end
                orPos = newPos(1, 1:2);
                set(dragging,'Position',current_position);
            else
                rectangle_position = get(a, 'Position');
                mouse_position = get(gca,'CurrentPoint');
                new_radius = (mouse_position(1) - rectangle_position(1));
                new_rectangle_position = [rectangle_position(1:2) new_radius new_radius];
                if all(new_rectangle_position > 0)
                    set(a, 'Position', new_rectangle_position)
                end
            end
        end
    end
  drawnow;
end


%%%%%% CALLBACKS %%%%%
function openDicom(hObject, eventdata)
    handles = guidata(hObject);
    
    if isfield(handles, 'lastPath')
        dirName = uigetdir(handles.lastPath, 'Select target DICOM folder');
    else
        dirName = uigetdir('Select target DICOM folder');
    end
    
    if dirName
        
        handles.lastPath = dirName;
        
        [metadata, dicomImages] = getDicomData(dirName);   
        
        lung = single(dicomImages);
        
        %Sort lung
        %[handles.data.lung, handles.data.metadata] =...
        %    sortLungOrMask(lung, metadata);
        
        handles.data.lung = lung;
        handles.data.metadata = metadata;
        
        %Post opening
        prepareAndShowLung(handles, 1)
        makeWidgetsVisible(handles)
        
        guidata(hObject, handles)
        
        %Set patient name
        set(handles.gui.patName, 'String', dirName);
    end
end

%Save Results
function saveResults(hObject, eventdata)
    handles = guidata(hObject);
    
    patientName = getPatientName(handles.data.metadata(1));
    
    [name, pathName] = uiputfile([patientName '.mat'],...
        'Save Results');
    
    if name
        %Whole Lung
        allResults.wholeLung.classical = handles.data.resultsClassicalWholeLung;
        allResults.wholeLung.percentual = handles.data.resultsPercentileWholeLung;
        %Base Lung
        allResults.baseLung.classical = handles.data.resultsClassicalBaseLung;      
        allResults.baseLung.percentual = handles.data.resultsPercentileBaseLung;
        %Middle Lung
        allResults.middleLung.classical = handles.data.resultsClassicalMiddleLung;      
        allResults.middleLung.percentual = handles.data.resultsPercentileMiddleLung;
        %Top Lung
        allResults.topLung.classical = handles.data.resultsClassicalTopLung;      
        allResults.topLung.percentual = handles.data.resultsPercentileTopLung;
        
        %Structure
        allResults.structure.mask = handles.data.lungMask;
        allResults.structure.lung = handles.data.lung;
        allResults.structure.roiAir = handles.data.roiAir;
        allResults.structure.roiTissue = handles.data.roiTissue;
        allResults.structure.metadata = handles.data.metadata;
        allResults.structure.errata = 'Changed poor to normally';
        save([pathName name], 'allResults')
    end
end

function openAirWay(hObject, eventdata)
    handles = guidata(hObject);
    
    [fileName pathName] = uigetfile('*.nrrd','Choose AIRWAY file',...
        handles.lastPath);
    
    if fileName
        masks = nrrd_read([pathName fileName]);
        airWayMask = flipdim(masks, 3);
        handles.data.lungMask(airWayMask  >= 1) = 0;
        guidata(hObject, handles)
    end
   
end

function refreshFunctionsList(hObject, eventdata)
    handles = guidata(hObject);
    setFunctionsName(handles);
end

function roiAirCallback(hObject, eventdata)
    handles = guidata(hObject);
    [handles.roi_air_properties.handle, handles.roi_air_properties.position] =...
        drag_and_drop(100, 100, 8, 'g'); 
    guidata(hObject, handles)
end

function roiTissueCallback(hObject, eventdata)
    handles = guidata(hObject);
    [handles.roi_tissue_properties.handle, handles.roi_tissue_properties.position] =...
        drag_and_drop(100, 100, 8, 'r'); 
    guidata(hObject, handles)
end

function openHDR(hObject, eventdata)
    handles = guidata(hObject);
    if isfield(handles, 'lastPath')
        [name pathName] = uigetfile('*.hdr',...
            'Select .HDR mask file',...
            handles.lastPath);
    else
        [name pathName] = uigetfile('*.hdr', 'Select .HDR mask file');
    end
    
    if name
        masks = getHDRMask([pathName name]);
        
                
        %Sort masks
        %handles.data.lungMask = sortLungOrMask(masks, handles.data.metadata);      
     
        handles.data.lungMask = masks;
        guidata(hObject, handles)
        
        %Enable widget On
        set(handles.gui.showMask, 'Visible', 'On')
        set(handles.gui.executeButton, 'Visible', 'On'); 
        %menu
        set(handles.gui.calibrationMenu, 'Enable', 'On')
        set(handles.gui.imageMenu, 'Enable', 'On')        
    end
end

function openNrrdMask(hObject, eventdata)
    handles = guidata(hObject);
    
    [fileName pathName] = uigetfile('*.nrrd','Choose AIRWAY file',...
        handles.lastPath);
    
    if fileName
        masks = nrrd_read([pathName fileName]);
        nTypes = unique(masks);
        
        %Chest Image Platform
        if length(nTypes) > 2
            masks(masks > 500) = 0;            
        end
        
        masks(masks >= 1) = 1;
        
        %Sort masks
        %handles.data.lungMask = sortLungOrMask(masks,...
        %   handles.data.metadata);    
        
        handles.data.lungMask = masks;
      
        set(handles.gui.showMask, 'Visible', 'On')
        set(handles.gui.imageMenu, 'Enable', 'On')
        set(handles.gui.calibrationMenu, 'Enable', 'On')
        set(handles.gui.executeButton, 'Visible', 'On')
        guidata(hObject, handles);
    end
end

function moveSlice(hObject, eventdata, direction)
    handles = guidata(hObject);
    lung = handles.data.lung;
    currentSlicePosition = getCurrentSlicePosition(handles);
    if direction
        currentSlicePosition = currentSlicePosition + 1;
    else
        currentSlicePosition = currentSlicePosition - 1;
    end
    
    newSlicePosition = decideNewSlice(lung, currentSlicePosition);
    setCurrentSlicePosition(handles, newSlicePosition)
    prepareAndShowLung(handles, newSlicePosition)
end

function newSlicePosition = decideNewSlice(lung, currentSlicePosition)
    nSlices = size(lung, 3);
    if currentSlicePosition > nSlices
        newSlicePosition = 1;
    elseif currentSlicePosition == 0
        newSlicePosition = nSlices;
    else
        newSlicePosition = currentSlicePosition;
    end
end

function flipImage(hObject, eventdata)
    handles = guidata(hObject);
    handles.data.lung = flipdim(handles.data.lung, 3);
    guidata(hObject, handles)
end

function flipMask(hObject, eventdata)
    handles = guidata(hObject);
    handles.data.lungMask = flipdim(handles.data.lungMask, 3);
    guidata(hObject, handles)
end


function execute(hObject, eventdata)
    handles = guidata(hObject);
    
    [hasROIS, avgAir, avgTissue, roiAir, roiTissue] =...
        checkForCalibrationROI(handles);
    
    if hasROIS
        handles.data.avgAir = avgAir;
        handles.data.avgTissue = avgTissue;
        

        parenchyma = handles.data.lung;
        mask = handles.data.lungMask;
        mask(mask >= 1) = 1;
        parenchyma(mask ~= 1) = 10000;
        
        
        if min(parenchyma(:)) < 0
            handles.data.lung = uncalibrateLung(handles.data.lung,...
                handles.data.metadata{1});
        end
        
        %Calibrate Lung
        handles.data.lung = lungCalibration(handles.data.lung,...
            handles.data.avgAir, handles.data.avgTissue);
        [handles.data.resultsClassicalWholeLung,...
            handles.data.resultsPercentileWholeLung,...
            handles.data.resultsClassicalBaseLung,...
            handles.data.resultsPercentileBaseLung,...
            handles.data.resultsClassicalMiddleLung,...
            handles.data.resultsPercentileMiddleLung,...
            handles.data.resultsClassicalTopLung,...
            handles.data.resultsPercentileTopLung]  =...
            allAnalysis(handles.data.lung, handles.data.lungMask, handles.data.metadata);
        
        handles.data.roiAir = roiAir;
        handles.data.roiTissue = roiTissue;
        handles.data.lung = handles.data.lung;
        %Enable save results menu
        set(handles.gui.saveResults, 'Enable', 'On');
        
        guidata(hObject, handles)
    end
    
end


%%%%%% INTERFACE %%%%%%
function drawInterface()
    screenSize = get(0, 'ScreenSize');
    mainFigure = figure('IntegerHandle','off',...
        'Menubar', 'None',...
        'Name', 'LIP: Lung Image Processing 0.0.3dev0',...
        'Position', screenSize,...
        'tag', 'mainFig');

    mainAxes = axes('Parent', mainFigure,...
       'tag', 'mainAxes',...
       'Visible', 'off');

    bckColor = get(mainFigure, 'Color');
    
    fileMenu = uimenu('Parent', mainFigure,...
        'Label', 'File');
    
    %% Menus %%
    uimenu('Parent', fileMenu,...
        'Label', 'Load DICOM',...
        'Callback', @openDicom);
    
    masksMenu = uimenu('Parent', fileMenu,...
        'Label', 'Open Masks',...
        'Tag', 'maskMenu',...
        'Enable', 'Off');
    
    uimenu('Parent', fileMenu,...
        'Label', 'Save Results',...
        'Tag', 'saveResults',...
        'Enable', 'Off',...
        'Callback', @saveResults);
    
    uimenu('Parent', masksMenu,...
        'Label', '.hdr',...
        'Callback', @openHDR)
    
   uimenu('Parent', masksMenu,...
        'Label', '.nrrd',...
        'Callback', @openNrrdMask);
    
%     uimenu('Parent', nrrdMenu,...
%         'Label', 'without airway');
%     
%     uimenu('Parent', nrrdMenu,...
%         'Label', 'with airway');
    
        
    uimenu('Parent', masksMenu,...
        'Label', 'Open airway',...
        'Callback', @openAirWay);
    
    calibrationMenu = uimenu('Parent',mainFigure,...
        'Label', 'Calibration',...
        'Tag', 'calibrationMenu',...
        'Enable', 'Off');
    
    calibROIMenu = uimenu('Parent', calibrationMenu,...
        'Label', 'Set Calibration ROIs');
    
    uimenu('Parent', calibROIMenu,...
        'Label', 'ROI Trachea',...
        'Callback', @roiAirCallback);
    
    uimenu('Parent', calibROIMenu,...
        'Label', 'ROI Aorta',...
        'Callback', @roiTissueCallback);
    
    imageMenu = uimenu('Parent', mainFigure,...
        'Label', 'Image',...
        'Tag', 'imageMenu',...
        'Enable', 'Off');
    
    uimenu('Parent', imageMenu,...
        'Label', 'Flip Image',...
        'Callback', @flipImage);    
        
    uimenu('Parent', imageMenu,...
        'Label', 'Flip Mask',...
        'Callback', @flipMask);
       
    
    %% Buttons %%
    uicontrol('Parent', mainFigure,...
        'Units', 'Normalized',...
        'Position',  [0.24, 0.1, 0.05, 0.03],...
        'String', 'Previous',...
        'Callback', {@moveSlice, 0},...
        'Tag', 'buttonPrevious',...
        'Visible', 'Off')
    
    uicontrol('Parent', mainFigure,...
        'Units', 'Normalized',...
        'Position',  [0.75, 0.1, 0.05, 0.03],...
        'String', 'Next',...
        'Callback', {@moveSlice, 1},...
        'Tag', 'buttonNext',...
        'Visible', 'Off')
    
    uicontrol('Parent', mainFigure,...
        'Units', 'Normalized',...
        'Position',[0.15, 0.1, 0.05, 0.03],...
        'String', 'Execute',...
        'Visible', 'Off',...
        'Tag', 'executeButton',...
        'Callback', @execute)    
    
   
    %% Edit box %%
    uicontrol('Parent', mainFigure,...
        'Units', 'Normalized',...
        'Position', [0.5, 0.07, 0.05, 0.03],...
        'Style', 'Edit',...
        'String', '1',...
        'Tag', 'currentSlicePosition',...
        'Backgroundcolor', [1 1 1],...
        'Visible', 'Off')
    
    
    %% Check box %%
    uicontrol('Parent', mainFigure,...
        'Units', 'Normalized',...
        'Position', [0.49, 0.025, 0.09, 0.03],...
        'Style', 'Check',...
        'String', 'Show Mask',...
        'Tag', 'showMask',...
        'Visible', 'Off')
    
    %% Static Text %%
    uicontrol('Parent', mainFigure,...
        'Units', 'Normalized',...
        'Position', [0.35, 0.93, 0.4, 0.05],...
        'Style', 'Text',...
        'HorizontalAlignment', 'Left',...
        'Backgroundcolor', bckColor,...
        'String', 'Patient Name:',...
        'Tag', 'patName',...
        'Visible', 'Off')    
        
    handles.gui = guihandles(mainFigure);
    guidata(mainFigure, handles)

end


%%% EXTERNAL FUNCTION %%%

function [X] = nrrd_read(filename)
% Modified from the below function (only X is outputed and cleaner is not used)
%NRRDREAD  Import NRRD imagery and metadata.
%   [X, META] = NRRDREAD(FILENAME) reads the image volume and associated
%   metadata from the NRRD-format file specified by FILENAME.
%
%   Current limitations/caveats:
%   * "Block" datatype is not supported.
%   * Only tested with "gzip" and "raw" file encodings.
%   * Very limited testing on actual files.
%   * I only spent a couple minutes reading the NRRD spec.
%
%   See the format specification online:
%   http://teem.sourceforge.net/nrrd/format.html

% Copyright 2012 The MathWorks, Inc.

try
    % Open file.
    fid = fopen(filename, 'rb');
    assert(fid > 0, 'Could not open file.');
%     cleaner = onCleanup(@() fclose(fid));

    % Magic line.
    theLine = fgetl(fid);
    assert(numel(theLine) >= 4, 'Bad signature in file.')
    assert(isequal(theLine(1:4), 'NRRD'), 'Bad signature in file.')

    % The general format of a NRRD file (with attached header) is:
    % 
    %     NRRD000X
    %     <field>: <desc>
    %     <field>: <desc>
    %     # <comment>
    %     ...
    %     <field>: <desc>
    %     <key>:=<value>
    %     <key>:=<value>
    %     <key>:=<value>
    %     # <comment>
    % 
    %     <data><data><data><data><data><data>...

    meta = struct([]);

    % Parse the file a line at a time.
    while (true)

      theLine = fgetl(fid);

      if (isempty(theLine) || feof(fid))
        % End of the header.
        break;
      end

      if (isequal(theLine(1), '#'))
          % Comment line.
          continue;
      end

      % "fieldname:= value" or "fieldname: value" or "fieldname:value"
      parsedLine = regexp(theLine, ':=?\s*', 'split','once');

      assert(numel(parsedLine) == 2, 'Parsing error')

      field = lower(parsedLine{1});
      value = parsedLine{2};

      field(isspace(field)) = '';
      meta(1).(field) = value;

    end

    datatype = getDatatype(meta.type);

    % Get the size of the data.
    assert(isfield(meta, 'sizes') && ...
           isfield(meta, 'dimension') && ...
           isfield(meta, 'encoding') && ...
           isfield(meta, 'endian'), ...
           'Missing required metadata fields.')

    dims = sscanf(meta.sizes, '%d');
    ndims = sscanf(meta.dimension, '%d');
    assert(numel(dims) == ndims);

    data = readData(fid, meta, datatype);
    data = adjustEndian(data, meta);

    % Reshape and get into MATLAB's order.
    X = reshape(data, dims');
    X = permute(X, [2 1 3]);
    
    fclose(fid);
catch err
    fclose(fid);
    rethrow(err)
end
end

function datatype = getDatatype(metaType)

% Determine the datatype
switch (metaType)
 case {'signed char', 'int8', 'int8_t'}
  datatype = 'int8';
  
 case {'uchar', 'unsigned char', 'uint8', 'uint8_t'}
  datatype = 'uint8';

 case {'short', 'short int', 'signed short', 'signed short int', ...
       'int16', 'int16_t'}
  datatype = 'int16';
  
 case {'ushort', 'unsigned short', 'unsigned short int', 'uint16', ...
       'uint16_t'}
  datatype = 'uint16';
  
 case {'int', 'signed int', 'int32', 'int32_t'}
  datatype = 'int32';
  
 case {'uint', 'unsigned int', 'uint32', 'uint32_t'}
  datatype = 'uint32';
  
 case {'longlong', 'long long', 'long long int', 'signed long long', ...
       'signed long long int', 'int64', 'int64_t'}
  datatype = 'int64';
  
 case {'ulonglong', 'unsigned long long', 'unsigned long long int', ...
       'uint64', 'uint64_t'}
  datatype = 'uint64';
  
 case {'float'}
  datatype = 'single';
  
 case {'double'}
  datatype = 'double';
  
 otherwise
  assert(false, 'Unknown datatype')
end
end


function data = readData(fidIn, meta, datatype)

switch (meta.encoding)
 case {'raw'}
  
  data = fread(fidIn, inf, [datatype '=>' datatype]);
  
 case {'gzip', 'gz'}

  tmpBase = tempname();
  tmpFile = [tmpBase '.gz'];
  fidTmp = fopen(tmpFile, 'wb');
  assert(fidTmp > 3, 'Could not open temporary file for GZIP decompression')
  
  tmp = fread(fidIn, inf, 'uint8=>uint8');
  fwrite(fidTmp, tmp, 'uint8');
  fclose(fidTmp);
  
  gunzip(tmpFile)
  
  fidTmp = fopen(tmpBase, 'rb');
%   cleaner = onCleanup(@() fclose(fidTmp));
  
  meta.encoding = 'raw';
  data = readData(fidTmp, meta, datatype);
  fclose(fidTmp)
  
 case {'txt', 'text', 'ascii'}
  
  data = fscanf(fidIn, '%f');
  data = cast(data, datatype);
  
 otherwise
  assert(false, 'Unsupported encoding')
end
end


function data = adjustEndian(data, meta)

[void,void,endian] = computer();

needToSwap = (isequal(endian, 'B') && isequal(lower(meta.endian), 'little')) || ...
             (isequal(endian, 'L') && isequal(lower(meta.endian), 'big'));
         
if (needToSwap)
    data = swapbytes(data);
end
end