function logSrepQuadPrim = Log(srepQuadPrim)

if( sum([srepQuadPrim.inTangentSpace]) ~= 0 )
	srepQuadPrim
	error( 'symmetric-space primitive is not in manifold space.' );
end

logSrepQuadPrim = srepQuadPrim;

for i=1:numel(srepQuadPrim)
    logSrepQuadPrim(i) = SrepQuadPrimitive(srepQuadPrim.pos, ...
        log(srepQuadPrim.r), SphereLog(srepQuadPrim.U));
end
