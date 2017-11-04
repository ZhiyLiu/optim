% ==============================================================
% XXX Need modification for multi-figure objects
% --------------------------------------------------------------
function ret = readLandmarks(m3dFilename)

    arrayChar = flat(m3dFilename);
    %nObjs = sscanf(arrayChar(strmatch('model.figureTrees.count', arrayChar), :),  'model.figureTrees.count %g');
    nObjs = findVal('model.figureTrees.count', arrayChar, ' %g');
    
    landmarks = cell(1, nObjs); % each landmark consists of 'name', [x y z], atomIndex. 
    for i=1:nObjs
        figKeyStr = [ 'model.figure[' int2str(i-1) '].'];
        keyStr = [figKeyStr 'numLandmarks'];
        %nMarks = sscanf(arrayChar(strmatch(keyStr, arrayChar), :), [keyStr ' %g']);
        nMarks = findVal(keyStr, arrayChar, ' %g');
        
        names = cell(1, nMarks);
        atomIndices = zeros(1, nMarks);
        points = zeros(3, nMarks);
        isPoints = zeros(1, nMarks);  %flags
        for j=1:nMarks
            markKeyStr = ['landmark[' int2str(j-1) '].'];
            keyStr = [figKeyStr markKeyStr 'name'];
            names{j} = findVal(keyStr, arrayChar, ' %s');            
            
            keyStr = [figKeyStr markKeyStr 'atomIndex'];       
            atomIndex = findVal(keyStr, arrayChar, ' %g');
            
            if (isempty(atomIndex))
                keyStr = [figKeyStr markKeyStr 'x'];       
                x = findVal(keyStr, arrayChar, ' %g');
                keyStr = [figKeyStr markKeyStr 'y'];       
                y  = findVal(keyStr, arrayChar, ' %g');
                keyStr = [figKeyStr markKeyStr 'z'];       
                z = findVal(keyStr, arrayChar, ' %g');
                if (~isempty(x) & ~isempty(y) & ~isempty(z))
                    points(:, j) = [x, y, z]';
                end
                isPoints(j) = 1;    %set true
            else
                atomIndices(j) = atomIndex;
            end
        end
        %numlandmarks, names, atomIndices, points
        landmarks{i} = {nMarks, names, atomIndices, points, isPoints};
    end
    
    ret = cell(1);
    ret{1} = landmarks;
return;
    
