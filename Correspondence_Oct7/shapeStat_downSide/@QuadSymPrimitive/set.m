function qSymPrim = set(qSymPrim,varargin)
% SET Set asset properties and return the updated object
propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);    
    switch prop
    case 'pos'
        qSymPrim.pos = val;
    case 'r'
        qSymPrim.r = val;
    case 'Up1'
        qSymPrim.Up1 = val;
    case 'Um1'
        qSymPrim.Um1 = val;        
    case 'elongation'
        qSymPrim.elongation = val;
    case 'inTangentSpace'
        qSymPrim.inTangentSpace = val;
    otherwise
        error([propName,' Is not a valid QuadSymPrimitive property'])
    end
end