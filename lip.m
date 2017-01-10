function lip()
    drawInterface()
end

function [metadata, dicomImages] = getDicomData(dirName)
    filesAndFolders = dir(dirName);
    nItems = length(filesAndFolders);
    
    
    h = waitbar(0,'Loading Metadata and DICOMS...');
    
    counter = 1;
    for idx = 1:nItems
        currentItemName = filesAndFolders(idx).name;
        if ~strcmp(currentItemName, '.') && ~strcmp(currentItemName, '..')...
                && ~strcmp(currentItemName, '.DS_Store')
            try
                fullPath = [dirName filesep currentItemName];
                metadata(counter) = dicominfo(fullPath);
                dicomImages(:, :, counter) = dicomread(fullPath);
                counter = counter + 1;
                waitbar(counter / nItems);
            catch
                continue
            end
        end
    end
    close(h)
end

%%%% LUNG CALIBRATION %%%%
function calibratedLung = lungCalibration(rawLung, mAir, mTissue)
    coef = polyfit([mAir, mTissue], [-1000 50], 1);
    calibratedLung = rawLung * coef(1) + coef(2);
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
    set(handles.gui.executeButton, 'Visible', 'On');
end

function lung = uncalibrateLung(lung, metadata)
    lung = (lung - metadata.RescaleIntercept) / metadata.RescaleSlope;
end


function [hasROIS, output] = checkForCalibrationROI(handles)
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
        output.avg.air = avgAir;
        output.avg.tissue = avgTissue;
        output.struct.air = roiAir;
        output.struct.tissue = roiTissue;
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
        handles.data.metadata = metadata;
        handles.data.lung = uncalibrateLung(dicomImages, metadata(1));
        guidata(hObject, handles)
        
        %Post opening
        prepareAndShowLung(handles, 1)
        makeWidgetsVisible(handles)
        
        %Set patient name
        set(handles.gui.patName, 'String', dirName);
    end
end

function openAirWay(hObject, eventdata)
    handles = guidata(hObject);
    
    [fileName pathName] = uigetfile('*.nrrd','Choose AIRWAY file',...
        handles.lastPath);
    
    if fileName
        masks = nrrd_read([pathName fileName]);
        airWayMask = flipdim(masks, 3);
        handles.data.lungMask(airWayMask >= 1) = 0;
         guidata(hObject, handles)
    end
   
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
        handles.data.lungMask = masks;
        guidata(hObject, handles)
        set(handles.gui.showMask, 'Visible', 'On')
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
    
    [hasROIS, roiInfo] = checkForCalibrationROI(handles);
    
    if hasROIS
        avgAir = roiInfo.avg.air;
        avgTissue = roiInfo.avg.tissue;
        
        %Calibrate Lung
        calibratedLung = lungCalibration(handles.data.lung, avgAir,...
            avgTissue);
    end
    
end


%%%%%% INTERFACE %%%%%%
function drawInterface()
    screenSize = get(0, 'ScreenSize');
    mainFigure = figure('IntegerHandle','off',...
        'Menubar', 'None',...
        'Name', 'Batch P15',...
        'Number', 'Off',...
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
        'Label', 'Open Masks');
    
    uimenu('Parent', masksMenu,...
        'Label', '.hdr',...
        'Callback', @openHDR)
    
    nrrdMenu = uimenu('Parent', masksMenu,...
        'Label', '.nrrd');
    
    uimenu('Parent', nrrdMenu,...
        'Label', 'without airway');
    
    uimenu('Parent', nrrdMenu,...
        'Label', 'with airway');
    
        
    uimenu('Parent', masksMenu,...
        'Label', 'Open airway',...
        'Callback', @openAirWay);
    
    calibrationMenu = uimenu('Parent',mainFigure,...
        'Label', 'Calibration');
    
    calibROIMenu = uimenu('Parent', calibrationMenu,...
        'Label', 'Set Calibration ROIs');
    
    uimenu('Parent', calibROIMenu,...
        'Label', 'ROI Trachea',...
        'Callback', @roiAirCallback);
    
    uimenu('Parent', calibROIMenu,...
        'Label', 'ROI Aorta',...
        'Callback', @roiTissueCallback);
    
    imageMenu = uimenu('Parent', mainFigure,...
        'Label', 'Image');
    
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