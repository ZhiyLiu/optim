function qP = set(qP,varargin)
% SET Set asset properties and return the updated object
propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);    
    switch prop
        case 'pos'
            qP.pos = val;
        case 'r'
            qP.r = val;
        case 'theta'
            qP.theta = val;
        case 'q'
            qP.q = val;
        case 'elongation'
            qP.elongation = val;
        case 'inTangentSpace'
            qP.inTangentSpace = val;
        otherwise
            error('QuadPrimitive properties: x, y, z, r, thetq, q, elongation, inTangentSpace');
    end
end
     
            
            
  
