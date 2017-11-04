% ==============================================================
% XXX Need modification for multi-figure objects
% --------------------------------------------------------------
function ret = readFigNames(m3dFilename)

    arrayChar = flat(m3dFilename);
    nObjs = sscanf(arrayChar(strmatch('model.figureTrees.count', arrayChar), :),  'model.figureTrees.count %g');
    
    figNames = cell(1, nObjs); % each landmark consists of 'name', [x y z], atomIndex. 
    for i=1:nObjs
        figKeyStr = [ 'model.figure[' int2str(i-1) '].'];
        keyStr = [figKeyStr 'name'];
        figNames{i} = findVal(keyStr, arrayChar, ' %s');
    end
    
    ret = cell(1);
    ret{1} = figNames; 
return;
