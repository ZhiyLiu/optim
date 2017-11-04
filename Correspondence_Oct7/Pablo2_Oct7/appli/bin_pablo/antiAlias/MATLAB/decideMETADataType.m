function matlabDataType = decideMETADataType(metaDataType)
%
% function matlabDataType = decideMETADataType(metaDataType)
%	A convience function to map META data type to matlab data type.
% Input
%	metaDataType	The string in the meta header file.
% Output
%	matlabDataType	The corresponding matlab data type.
%
if strcmpi('MET_FLOAT',metaDataType)
	matlabDataType = 'float';
elseif strcmpi('MET_DOUBLE',metaDataType)
	matlabDataType = 'float64';
elseif strcmpi('MET_USHORT',metaDataType)
	matlabDataType = 'uint16';
elseif strcmpi('MET_SHORT',metaDataType)
	matlabDataType = 'int16';
elseif strcmpi('MET_UCHAR',metaDataType)
	matlabDataType = 'uint8';
else
	error(['Unknown data type <' metaDataType '>']);
end