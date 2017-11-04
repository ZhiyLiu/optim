function antiAliasWrapper( inFile, outPrefix, movement )
%
% Anti-aliases. This function is just a wrapper around the
% antiAliasFunction being used.
%
% function antiAliasWrapper( inFile, outPrefix, movement )
%	Anti aliases a 3D binary image (1 inside, 0 outside and the boundary
% assumed to lie at 0.5); return a level-set with 0 representing the
% boundary. Movement should be [0.5 0.5] for strict antiAlias, anything
% more will result in smoothing. The first value is for the internal
% threshold, the second is for the outer threshold.
%

I = readMetaImage( inFile );

%[antiAliased, K, H] = antiAlias3D_dH2( I.I, 1, movement, 4, 2000, I.spacing );

% xiaojie--use large bandwidth -0425
bandwidth = 6; %min(max(ceil(movement*2)),3);

[antiAliased, K, H] = antiAlias3D_dH_flatRegions( I.I, 1, movement, bandwidth, 2000, I.spacing );
%[antiAliased, K, H] = antiAlias3D_dH( I.I, 1, movement, 2000, I.spacing );

% invert the signed distance map before saving so that marching cubes
% doesn't produce an inverted tile set.

% xiaojie ---shift then time 100 for antiAliased 0421
I.I = antiAliased*100 + 65536*max(-sign(antiAliased),0);
% I.I =(antiAliased)*100;
% I.I =(antiAliased-min(min(min(antiAliased))))*100;%xiaojie 04/16/11-- consistant with Pablo ImageDistanceMap.cpp
%xiaojie 04/16/11--add'MET_USHORT',
% get the results used by mhd2raw3.pl
writeMetaImage( I, [outPrefix '-antiAliased.mhd'],'MET_USHORT');

% [I1, I2, I3] = gradient(I.I);
% I.I = I1;
% writeMetaImage( I, [outPrefix '-antiAliased_gradX.mhd'],'MET_USHORT' );
% I.I = I2;
% writeMetaImage( I, [outPrefix '-antiAliased_gradY.mhd'],'MET_USHORT' );
% I.I = I3;
% writeMetaImage( I, [outPrefix '-antiAliased_gradZ.mhd'],'MET_USHORT' );
% 
% I.I = K;
% writeMetaImage( I, [outPrefix '-K.mhd'],'MET_USHORT' );
% I.I = H;
% writeMetaImage( I, [outPrefix '-H.mhd'],'MET_USHORT');

end
