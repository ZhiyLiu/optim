function val = get(srepQuadPrim, propName)

switch propName
    case 'pos'
        val = [srepQuadPrim.pos];
    case 'r'
        val = [srepQuadPrim.r];
    case 'U'
        val = [srepQuadPrim.U];
    case 'inTangentSpace'
        val = [srepQuadPrim.inTangentSpace];
    case 'Type'
        val = 's-rep' ;        
    otherwise
        error([propName,' Is not a valid SrepQuadPrimitive property'])
end
