function  simTran = parseSimTrans(modelArray)
    
    scale = findVal('model.transformation.scale', modelArray, ' %f');
    rw = findVal('model.transformation.rotation.w', modelArray, ' %f');
    rx = findVal('model.transformation.rotation.x', modelArray, ' %f');
    ry = findVal('model.transformation.rotation.y', modelArray, ' %f');
    rz = findVal('model.transformation.rotation.z', modelArray, ' %f');
    tx = findVal('model.transformation.translation.x', modelArray, ' %f');
    ty = findVal('model.transformation.translation.y', modelArray, ' %f');
    tz = findVal('model.transformation.translation.z', modelArray, ' %f');
	
	% FIXME: Need to include new simtrans info in here?

    simTran = [tx; ty; tz; rw; rx; ry; rz; scale];    
    
return;
