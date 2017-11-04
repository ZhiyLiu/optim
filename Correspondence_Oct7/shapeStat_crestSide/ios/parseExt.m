function worldExt = parseExt(modelArray)
    baseStr = 'model.world.';
    imageModTime = findVal([baseStr 'imageModTime'], modelArray, ' %s %s %s %s %s'); 
    imagePath = findVal([baseStr 'imagePath'], modelArray, ' s'); 
    
    bX = findVal([baseStr 'bound.x'], modelArray, ' %f'); 
    if (isempty(bX))
       worldExt = {};
       return;
    end
    
    bY = findVal([baseStr 'bound.y'], modelArray, ' %f'); 
    bZ = findVal([baseStr 'bound.z'], modelArray, ' %f'); 
    oX = findVal([baseStr 'origin.x'], modelArray, ' %f'); 
    oY = findVal([baseStr 'origin.y'], modelArray, ' %f'); 
    oZ = findVal([baseStr 'origin.z'], modelArray, ' %f'); 
    sX = findVal([baseStr 'spacing.x'], modelArray, ' %f'); 
    sY = findVal([baseStr 'spacing.y'], modelArray, ' %f'); 
    sZ = findVal([baseStr 'spacing.z'], modelArray, ' %f'); 
    eX = bX - oX; 
    eY = bY - oY; 
    eZ = bZ - oZ;
    worldExt = {imageModTime, imagePath, [bX, bY, bZ], [oX, oY, oZ], [sX, sY, sZ], [eX, eY, eZ]};
return;
