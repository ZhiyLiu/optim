#ifndef _SURFACE_COLOR_MAP_H_
#define _SURFACE_COLOR_MAP_H_



/*
	This class stores or generates a colormap (greyscale or color ramp) for use
	by the Optimization Visualizer, class M3DObjectSurfaceVisualizer.

    This class supports several kinds of maps, which can be selected using
	function specifyMap().  The default is a linear greyscale map.

	A single colormap may be read from a file using the loadFromFile() function,
	in which case each line of the file should have a single R, G, and B triple.
    Colormaps are expected to be arrays with values ranging from 0 to 1 for red,
	green and blue, respectively.

	A capability for switching rendering among the loaded colormap and the
	generated maps is provided.
*/


class SurfaceColorMap
{
	public:

		// The names in P3DUserInterfaceCallback.cpp must correspond to these
		enum color_map_type {
			LoadedMap, GreyMap, InverseGreyMap, MagentaGreenMap
		};

		SurfaceColorMap();
		~SurfaceColorMap() { }

		// Load a colormap from a file.  This automatically sets
		// the loaded map to be the current map.
		void loadFromFile(const char * path);

		// Select the current colormap.  When this is called, the value
		// of P3DUserInterface::matchSurfaceColormapChoice should be adjusted.
		void specifyMap(color_map_type type);

		color_map_type currentMap() const { return t; }

//	protected:

		// Functions to access the current map.  These are used
		// in classes M3DObjectSurfaceVisualizer and P3View.

		// Length of the current colormap
		static int size() { return length; }

		// Component arrays of the current colormap
		std::vector<float> * red();
		std::vector<float> * green();
		std::vector<float> * blue();

	private:

		static int length;
		static std::vector<float> redVec;
		static std::vector<float> greenVec;
		static std::vector<float> blueVec;
		static std::vector<float> loadedRed;
		static std::vector<float> loadedGreen;
		static std::vector<float> loadedBlue;
		static color_map_type t;

		void initMagentaToGreen();
		void initGrey();
		void initInverseGrey();
};


// The only instance of class SurfaceColorMap
extern SurfaceColorMap surfaceColorMap;



#endif	/* _SURFACE_COLOR_MAP_H_ */

