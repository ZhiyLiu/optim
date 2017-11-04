function srepQuadPrim = set(srepQuadPrim, varargin)
% SET Set asset properties and return the updated object

propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);    
    switch prop
    case 'pos'
        srepQuadPrim.pos = val;
    case 'r'
        srepQuadPrim.r = val;
    case 'U'
        srepQuadPrim.U = val;
    case 'inTangentSpace'
        srepQuadPrim.inTangentSpace = val;
    otherwise
        error([propName,' Is not a valid QuadSymPrimitive property'])
    end
end
