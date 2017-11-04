function [parentFigIds rootFigIds] = readTree(modelArrays, nObjs, nChilds)


parentFigIds = zeros(1,numel(nChilds));
rootFigIds = zeros(1,numel(nChilds));
for i = 1:nObjs
    
    treeKeyStr = ['model.figureTrees.tree[' num2str(i-1) '].'];
    rootFigId = findVal([treeKeyStr 'figureId'], modelArrays, ' %d');
    depth = 0;
      
    [parentFigIds rootFigIds] = findChildFigId(modelArrays, treeKeyStr, nChilds, rootFigId, parentFigIds, rootFigIds, depth);
  
end    

    function  [parentFigIds rootFigIds] = findChildFigId(modelArrays, parentKeyStr, nChilds, parentFigId, parentFigIds, rootFigIds, depth)

        nSubFigs = nChilds(parentFigId+1);
        if (nSubFigs == 0)
            parentFigIds(parentFigId+1) = parentFigId;
            if (depth == 0)
                rootFigIds(parentFigId+1) = parentFigId;
            end
            return;
        end
        
        for k = 1:nSubFigs
            childKeyStr = [parentKeyStr 'child[' num2str(k-1) '].'];
            childFigId = findVal([childKeyStr 'figureId'], modelArrays, ' %d');
            
            idx = childFigId + 1;            
            parentFigIds(idx) = parentFigId;
            
            if (nChilds(idx) > 0)
                 depth = depth + 1;
                [parentFigIds rootFigIds]= findChildFigId(modelArrays, childKeyStr, nChilds, childFigId, parentFigIds, rootFigIds,depth);
            end

            if (depth == 0)
                rootFigIds(idx) = parentFigId;
            end
        end
    end

end