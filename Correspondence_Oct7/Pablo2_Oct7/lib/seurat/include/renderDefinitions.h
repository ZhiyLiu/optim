#ifndef RENDER_DEFINITIONS_H
#define RENDER_DEFINITIONS_H


// Value passed to Pointlist_server2::init().
extern const int MEDIAL_SUBDIVISIONS;			/* Plain tiles */
extern const int SUBDIV_MEDIAL_SUBDIVISIONS;	/* Subdivision tiles */

// Value passed to most calls of
// Pointlist_server2::ComputeTilePointCloud().
extern const int POINT_CLOUD_SUBDIVISIONS;

// Value passed to Pointlist_server2::ComputeTilePointCloud()
// in class Match.
extern const int MATCH_POINT_CLOUD_SUBDIVISIONS;

// Multiplier of level-of-detail value in surface renderer.  This
// times the level is passed to Pointlist_server2::init().
extern const int SURFACE_RESOLUTION_FACTOR;

// Offset of level-of-detail value in constraints rendering.  This plus
// the level is passed to Pointlist_server2::ComputeSubdivPointCloud().
extern const int CONSTRAINTS_RESOLUTION_BIAS;

// Used in M3DObjectRenderer.cpp.
// Must correspond to slider maxima in P3DUserInterface.fl.
extern const int ATOM_VECTORS_WIDTH_MAX;
extern const int MESH_CONNECTORS_WIDTH_MAX;

#endif

