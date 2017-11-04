#include <fstream>
#include <iostream>
#include <vector>

#include "SurfaceColorMap.h"


int SurfaceColorMap::length;
SurfaceColorMap::color_map_type SurfaceColorMap::t;

std::vector<float> SurfaceColorMap::redVec;
std::vector<float> SurfaceColorMap::greenVec;
std::vector<float> SurfaceColorMap::blueVec;

std::vector<float> SurfaceColorMap::loadedRed;
std::vector<float> SurfaceColorMap::loadedGreen;
std::vector<float> SurfaceColorMap::loadedBlue;

// Global instance of this class.  Note that this
// must be defined after the class variables above.
SurfaceColorMap surfaceColorMap;


std::vector<float> * SurfaceColorMap::red()
{
	if (t == LoadedMap)
		return &loadedRed;
	else
		return &redVec;
}

std::vector<float> * SurfaceColorMap::green()
{
	if (t == LoadedMap)
		return &loadedGreen;
	else
		return &greenVec;
}

std::vector<float> * SurfaceColorMap::blue()
{
	if (t == LoadedMap)
		return &loadedBlue;
	else
		return &blueVec;
}

SurfaceColorMap::SurfaceColorMap()
{
	specifyMap(GreyMap);	// Default
}

void SurfaceColorMap::specifyMap(color_map_type type)
{
	t = type;
	switch (type) {
		case LoadedMap:
			length = loadedRed.size();
			return;
		case GreyMap:
			initGrey();
			return;
		case InverseGreyMap:
			initInverseGrey();
			return;
		case MagentaGreenMap:
			initMagentaToGreen();
			return;
	}
}

void SurfaceColorMap::loadFromFile(const char * path)
{
	std::ifstream in(path);
	if (in) {
		loadedRed.clear();
		loadedGreen.clear();
		loadedBlue.clear();
		float r, g, b;
		while (in) {
			in >> r >> g >> b;
			if (in) {
				loadedRed.push_back(r);
				loadedGreen.push_back(g);
				loadedBlue.push_back(b);
				in.ignore(1000, '\n');
			}
		}
	}
	length = loadedRed.size();
	t = LoadedMap;
}

// [1, 0, 1] -> [0, 1, 0]
void SurfaceColorMap::initMagentaToGreen()
{
	redVec.clear();
	greenVec.clear();
    blueVec.clear();
	for (float i = 0; i <= 1.0f; i += 0.01f) {
		redVec.push_back(1.0f - i);
		greenVec.push_back(i);
		blueVec.push_back(1.0f - i);
	}
	length = redVec.size();
}

// [0, 0, 0] -> [1, 1, 1]
void SurfaceColorMap::initGrey()
{
	redVec.clear();
	greenVec.clear();
    blueVec.clear();

	redVec.push_back(0.0f);
	greenVec.push_back(0.0f);
	blueVec.push_back(0.0f);

	redVec.push_back(1.0f);
	greenVec.push_back(1.0f);
	blueVec.push_back(1.0f);

	length = redVec.size();
}

// [1, 1, 1] -> [0, 0, 0]
void SurfaceColorMap::initInverseGrey()
{
	redVec.clear();
	greenVec.clear();
    blueVec.clear();

	redVec.push_back(1.0f);
	greenVec.push_back(1.0f);
	blueVec.push_back(1.0f);

	redVec.push_back(0.0f);
	greenVec.push_back(0.0f);
	blueVec.push_back(0.0f);

	length = redVec.size();
}


