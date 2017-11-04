#include "M3DFigureStats.h"



/*
    This class is simply a container for statistics gleaned from the .m3d
    file, for use later in building the mask and computing image match.
    See objectfile/readFigureStats for how the data gets into this structure.
*/

M3DFigureStats::~M3DFigureStats()
{
	templates.clear();
	template_types.clear();
	rmses.clear();
	mean_offset.clear();
	mean_std.clear();
}

void M3DFigureStats::setTemplate(int template_ind, int sample, double val)
//sets the a sample of a template to val.
{
	int index;

    index = template_ind * dimension + sample;
    if(index < 0 || index > templates.size())
        return;

	if (index == templates.size())
		templates.push_back(val);
	else templates[index] = val;
}

void M3DFigureStats::setTemplateType(int point_ind, int template_ind)
//set the index to templates for the point point_ind.
{
	if (point_ind > template_types.size() || point_ind < 0)
		return;

	if (template_ind < 0 || template_ind >= numTemplates)
		return;

	if (point_ind == template_types.size())
		template_types.push_back(template_ind);
	else template_types[point_ind] = template_ind;
}

void M3DFigureStats::setRmses(int point_ind, double val)
//set the trained rms used to normalize point_ind profile during
//optimization.
{
	if (point_ind > rmses.size() || point_ind < 0)
		return;

	if (point_ind == rmses.size())
		rmses.push_back(val);
	else rmses[point_ind] = val;
}

void M3DFigureStats::setMeanOffset(int point_ind, double val)
//set mean intensity at point point_ind.  used for mean offset
//image match penalty.
{
	if (point_ind > mean_offset.size() || point_ind < 0)
		return;

	if (point_ind == mean_offset.size())
		mean_offset.push_back(val);
	else mean_offset[point_ind] = val;
}

void M3DFigureStats::setMeanStd(int point_ind, double val)
//set mean intensity standard deviation at point point_ind.  
//Used for mean offset image match penalty.
{
	if (point_ind > mean_std.size() || point_ind < 0)
		return;

	if (point_ind == mean_std.size())
		mean_std.push_back(val);
	else mean_std[point_ind] = val;
}

void M3DFigureStats::getTemplate(int point_ind, double v[])
//record the template required for point_ind point in the 
//array argument.
{
	int i, st = template_types[point_ind]*dimension;

	for (i = st; i < st + dimension; i ++) 
	{
		v[i] = templates[i];
	}
}

double M3DFigureStats::getTemplateVal(int point_ind, int sample)
//return only the element of the template for the point_ind
//point.
{
	int index;

    index = template_types[point_ind] * dimension + sample;
    if(index < 0 || index > templates.size())
        return -1.0;
	else return templates[index];
}

double M3DFigureStats::getTemplateValbyType(int temp_type, int sample)
//return template value given by template type and sample, rather
//than by point.
{
	int index;

	index = temp_type * dimension + sample;
	if (index < 0 || index > templates.size())
		return -1.0;
	else return templates[index];
}


