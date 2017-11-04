function quadPrims = convert2Lie(qSrepPrims)
% quadSimPrims = convert2Lie(srepQuadPrims)
%
% Convert Mx-format (s-rep) primitives to Lie-group representation
% primitives in the UNC format.  This is a lossy conversion since the
% opposing spokes may have different radii in the s-rep representation 
% but are forced to have the same radius in the UNC format.  The geometric 
% mean is used if they are different in the input.  
%
% The script is much shorter than this comment.

quadPrims = convert2Lie(convert2Sym(qSrepPrims));
