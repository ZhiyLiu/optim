function writeMetaImage(img, filename, dataType, extraTags)
%
% function [img] = writeMetaImage(img, filename, dataType, extraTags)
%	Writes 3D meta images and vector fields.
% Input
%	img			The img in the same structure form as returned by
%				readMetaImage.
%	filename	The file name of the meta image.
%	dataType	The meta header data type.
%	extraTags	A cell array of any additional tags that you want to pass.
% Output
%	none
%
if nargin < 4
	extraTags = {};
end
if nargin < 3
	dataType = 'MET_FLOAT';
end

[path,file,ext]=fileparts(filename);
%fprintf('path: %s\n', path);
%fprintf('file: %s\n', file);
%fprintf('ext: %s\n', ext);
if strcmp(ext,'') || strcmp(ext,'.mhd')
	if ~isempty(path)
		headerFilename = [path filesep file '.mhd'];
		dataFilename = [path filesep file '.raw'];
		dataFilenameHeaderRelative = [file '.raw'];
	else
		headerFilename = [file '.mhd'];
		dataFilename = [file '.raw'];
		dataFilenameHeaderRelative = [file '.raw'];
	end
else
	headerFilename = [filename '.mhd'];
	dataFilename = [filename '.raw'];
	dataFilenameHeaderRelative = [path filesep file '.' ext '.raw'];
end

%fprintf('header filename: %s\n', headerFilename);
%fprintf('data filename: %s\n', dataFilename);
%fprintf('data relative filename: %s\n', dataFilenameHeaderRelative);

nImDims = ndims(img.I);
numChannels = img.channels;
imSize = size(img.I);
field  = numChannels > 1;
if(field)
	nImDims = nImDims - 1;
	imSize  = imSize(1:end-1);
end

% write header file
%fprintf('writing header file: %s...\n',headerFilename);
fid = fopen(headerFilename,'w');
fprintf(fid,'ObjectType = %s\n','Image');
if field
	fprintf(fid,'ObjectSubType = %s\n','HField');
end
fprintf(fid,'NDims = %d\n',nImDims);
fprintf(fid,'BinaryData = %s\n','True');
fprintf(fid,'BinaryDataByteOrderMSB = %s\n','False');
fprintf(fid,'Offset = '); fprintf(fid,'%g ',img.origin); fprintf(fid,'\n');
fprintf(fid,'ElementSpacing = '); fprintf(fid,'%g ',img.spacing); fprintf(fid,'\n');
fprintf(fid,'DimSize = '); fprintf(fid,'%g ',imSize); fprintf(fid,'\n');
fprintf(fid,'ElementNumberOfChannels = %d\n',numChannels);
fprintf(fid,'ElementType = %s\n',dataType);
fprintf(fid,'ElementDataFile = %s\n',dataFilenameHeaderRelative);
for tagidx = 1:2:length(extraTags)
	fprintf(fid,'%s = %s\n', extraTags{tagidx}, extraTags{tagidx+1});
end

fclose(fid);

fid = fopen(dataFilename,'w');
%fprintf('writing data file: %s...\n',dataFilename);
if field
	% we need to permute the image before writing so that the channel
	% information is all together.
	I = permute(img.I, [nImDims+1 1:nImDims]);
else
	I = img.I;
end
count = fwrite(fid,I,decideMETADataType(dataType));
fclose(fid);
