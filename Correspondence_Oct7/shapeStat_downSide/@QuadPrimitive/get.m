function val = get(qPrim, propName)

switch propName
    case 'pos'
        val = qPrim.pos;
    case 'r'
        val = qPrim.r;
    case 'q'
        val = qPrim.q;
    case 'elongation'
        val = qPrim.elongation;
    case 'theta'
        val = qPrim.theta;
    case 'inTangentSpace'
        val = qPrim.inTangentSpace;
    otherwise
        error([propName,' Is not a valid QuadPrimitive property'])
end
