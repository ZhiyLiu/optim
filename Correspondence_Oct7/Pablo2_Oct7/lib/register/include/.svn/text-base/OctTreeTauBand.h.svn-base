/**
 * This is a data structure for storing (X, Y, Z) coordinates in
 * an oct-tree for easy search.
 *
 * This particular oct-tree is tailored for tau-band distance calculations
 * so it actually stores "original" and "normalized" coordinates
 * 
 * tau-band is defined in Josh Levy's dissertation as the Euclidean distance
 * between points that have acceptably close normalized coordinates 
 * and infinitely large otherwise.   It is used to enforce loose 
 * correspondence between the boundary of an m-rep and a set of 
 * manually specified points
 *
 * (c) 2008 Morphormics, Inc.
 */

#ifndef OCT_TREE_TAU_BAND_H
#define OCT_TREE_TAU_BAND_H

#include "Vector3D.h"


class OctTreeTauBand {

public:
  /**
   * A new tree needs an initial point (non-normalized)
   */
  OctTreeTauBand(const Vector3D& point);
  ~OctTreeTauBand();

  /**
   * Adding a node to the tree
   */
  OctTreeTauBand* insert (const Vector3D& point);

  /**
   * Calculate normalized coordinates relative to the extent of all the
   * points in the tree.  Note: if this function is never called, 
   * tau-band distance to that normalized to (0,0,0)  is really 
   * the Euclidean distance
   */
  void normalize();
  void normalize(const Vector3D& minExtent, const Vector3D& lengths);

  /**
   * Calculate the tau-band distance between a point and this tree.
   * If normalize() hasn't yet been called, this is really finding
   * the minimum squared Euclidean distance from the point to this tree.
   *
   * tau-band distance to points for which all normalized coordinates
   * are within tau of each other is Euclidean distance.
   *
   * tau-band distance to any other point is infinity
   *
   */
  double getDistanceSqr(const Vector3D& point, 
                        const Vector3D& normalizedPoint, double tau);

  /**
   * Calculate the tau-band distance from another tree to this tree.
   * If both trees have not yet been normalized, This is really finding
   * the sum of minimum squared Euclidean distances from the points in
   * the other tree to the this tree.  If only one tree has been normalized,
   * the result are unpredictable.
   */
  double getDistanceSqr(const OctTreeTauBand* other, double tau);

  /**
   * Calculate the average tau-band distance from another tree to this tree.
   * If both trees have not yet been normalized, This is really finding
   * the sum of minimum squared Euclidean distances from the points in
   * the other tree to the this tree.  If only one tree has been normalized,
   * the result are unpredictable.
   */
  double getAvgDistanceSqr(const OctTreeTauBand* other, double tau);

  /**
   * Debugging flags to control which distances the tau-band is enforced on
   */
  bool xLimit, yLimit, zLimit;

protected:

  /**
   * Defining the 8 possible neighbors for a voxel
   * Axial (z) coordinates are on the Inferior/Superior axis
   * Coronal (y) coordinates are on the Anterior/Posterior axis
   * Sagittal (x) coordinates are on the Left/Right axis
   */
  enum NeighborNames { RightAntInf, RightAntSup, RightPostInf, RightPostSup,
                         LeftAntInf, LeftAntSup, LeftPostInf, LeftPostSup };

  enum NeighborFlags { Sup = 1, Post = 2, Left = 4 };

  OctTreeTauBand* children[8];

  Vector3D pt; // in original coordinates

  Vector3D npt; // in normalized coordinates;

  /**
   * Keep track of the minimum and maximum extents of the tree
   * These values are only correct if insert() is never called
   * on a child without first being called on a parent
   */
  double minX, minY, minZ, maxX, maxY, maxZ;

  /**
   * Keep track of the number of nodes in the tree
   * This value is only correct if insert() is never called
   * on a child without first being called on a parent
   */
  int size;

  /**
   * really implements getDistanceSqr
   *
   * distance is the best distance seen so far
   */
  void getDistanceSqr(const Vector3D& point, 
                        const Vector3D& normalizedPoint, double tau,
			double& distance);

};

#endif /* OCT_TREE_H */


