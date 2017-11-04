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

#include "OctTreeTauBand.h"
#include <iostream>

/**
 * A new tree needs an initial point (non-normalized)
 */
OctTreeTauBand::OctTreeTauBand(const Vector3D& point) {

  // set the point coordinates
  pt.set(point);

  // every node starts out with the same normalized coordinates
  // so if you compute distances without normalizing first
  // it will give you Euclidean distances
  npt.set(0.0, 0.0, 0.0);

  // set the extents to be the only point in the tree
  minX = maxX = point.getX();
  minY = maxY = point.getY();
  minZ = maxZ = point.getZ();

  // I don't have any children yet;
  for (int child = 0; child < 8; child++) {
     children[child] = 0;
  }

  size = 1;

  xLimit = yLimit = false;
  zLimit = true;

}

/**
 * Cleanup... delete any children we've created
 */
OctTreeTauBand::~OctTreeTauBand() {
  for (int child = 0; child < 8; child++) {
    if (children[child]) {
      delete children[child];
      children[child] = 0;
    }
  }
}


/**
 * Adding a node to the tree
 */
OctTreeTauBand* OctTreeTauBand::insert (const Vector3D& point) {
  int whichNeighbor = 0;
  double x = point.getX();
  double y = point.getY();
  double z = point.getZ();

  // Figure out which neighbor this point wants to be
  // and update the max extents as we go along
  if (x < pt.getX()) {
    if (x < minX) {
      minX = x;
    }
    whichNeighbor |= Left;
  }
  else { // x >= pt.getX
    if (x > maxX) {
      maxX = x;
    }
  }

  if (y < pt.getY()) {
    if (y < minY) {
      minY = y;
    }
    whichNeighbor |= Post;
  }
  else { // y >= pt.getY
    if (y > maxY) {
      maxY = y;
    }
  }

  if (z < pt.getZ()) {
    if (z < minZ) {
      minZ = z;
    }
    whichNeighbor |= Sup;
  }
  else { // z >= pt.getZ
    if (z > maxZ) {
      maxZ = z;
    }
  }

  size++;

  // Once we've decided which child the voxel belongs to
  // Test to see if that child exists.  If so recurse on
  // insert.  If not create a new child.
  if (children[whichNeighbor]) {
    return children[whichNeighbor]->insert(point);
  }
  else {
    children[whichNeighbor] = new OctTreeTauBand(point);
	return children[whichNeighbor];
  }


}


/**
 * Calculate normalized coordinates relative to the extent of all the
 * points in the tree.  Note: if this function is never called, 
 * tau-band distance to that normalized to (0,0,0)  is really 
 * the Euclidean distance
 */
void OctTreeTauBand::normalize() {
  Vector3D minExtent(minX, minY, minZ);
  Vector3D lengths(maxX - minX, maxY - minY, maxZ - minZ);
  normalize(minExtent, lengths);
}

void OctTreeTauBand::normalize(const Vector3D& minExtent, const Vector3D& lengths) {
  npt.set(pt);
  npt -= minExtent;
  npt.setX(npt.getX() / lengths.getX());
  npt.setY(npt.getY() / lengths.getY());
  npt.setZ(npt.getZ() / lengths.getZ());

  for (int child = 0; child < 8; child++) {
    if (children[child]) {
      children[child]->normalize(minExtent, lengths);
    }
  }
}


/**
 * Calculate the tau-band distance between a point and this tree.
 * If normalize() hasn't yet been called, this is really finding
 * the minimum squared Euclidean distance from the point to this tree.
 *
 * tau-band distance to points for which all normalized coordinates
 * are within tau of each other is Euclidean distance.
 *
 * tau-band distance to any other point is infinity
 */
double OctTreeTauBand::getDistanceSqr(const Vector3D& point, 
                      const Vector3D& normalizedPoint, double tau) {
  
  double distance = 1000; //something huge compared to a unit cube;
  getDistanceSqr(point, normalizedPoint, tau, distance);
  return distance;
}

/**
 * Calculate the tau-band distance from another tree to this tree.
 * If both trees have not yet been normalized, This is really finding
 * the sum of minimum squared Euclidean distances from the points in
 * the other tree to the this tree.  If only one tree has been normalized,
 * the result are unpredictable.
 */
double OctTreeTauBand::getDistanceSqr(const OctTreeTauBand* other, double tau) {
  if (other) {
    double distSqr = getDistanceSqr(other->pt, other->npt, tau);
    for (int child = 0; child < 8; child++) {
      if (other->children[child]) {
        distSqr += getDistanceSqr(other->children[child], tau);
      }
    }
    return distSqr;
  }
  else {
    // distance from an empty set is 0
    return 0.0;
  }
}

/**
 * Calculate the average tau-band distance from another tree to this tree.
 * If both trees have not yet been normalized, This is really finding
 * the sum of minimum squared Euclidean distances from the points in
 * the other tree to the this tree.  If only one tree has been normalized,
 * the result are unpredictable.
 */
double OctTreeTauBand::getAvgDistanceSqr(const OctTreeTauBand* other, double tau) {
  if (other) {
    double sum = getDistanceSqr(other, tau);
    return sum / other->size;
  }
  else {
    return 0.0;
  }
}


/**
 * really implements getDistanceSqr
 *
 * distance is the best distance seen so far
 */
void OctTreeTauBand::getDistanceSqr(const Vector3D& point, 
                        const Vector3D& normalizedPoint, double tau,
			double& distance) {

  bool possibilities[8];
  for (int i=0; i<8; i++) {
    possibilities[i] = true;
  }

  Vector3D ndiff = normalizedPoint - npt;
  double ndx = ndiff.getX();
  double ndy = ndiff.getY();
  double ndz = ndiff.getZ();

  // double x = pt.getX();
  // double y = pt.getY();
  // double z = pt.getZ();

  // double nx = normalizedPoint.getX();
  // double ny = normalizedPoint.getY();
  // double nz = normalizedPoint.getZ();

  double tauSqr = tau*tau;

  double dist = pt.distSquare(point);

//  if (xLimit) {
//    std::cout << "xLimit is true" << std::endl;
//  }

  if (
      (xLimit && ((ndx * ndx) > tauSqr)) ||
      (yLimit && ((ndy * ndy) > tauSqr)) ||
      (zLimit && ((ndz * ndz) > tauSqr))
     ) {
     //dist += 1;
     //dist *= 1000;
     dist = 1000;
  }


  if (dist < distance) {
    distance = dist;
  }

  // use the oct tree wisely to search the subspace
  
  if (xLimit) {
    if (ndx > tau) { // the point I'm looking for is far to the right of where
                     // I'm at, so don't allow a search to the left
      possibilities[LeftAntInf] = false;
      possibilities[LeftAntSup] = false;
      possibilities[LeftPostInf] = false;
      possibilities[LeftPostSup] = false;
    }
    else if (ndx < -tau) { // the point I'm looking for is far to the left
                           // so disallow search to the right
      possibilities[RightAntInf] = false;
      possibilities[RightAntSup] = false;
      possibilities[RightPostInf] = false;
      possibilities[RightPostSup] = false;
    }
  }

  // Now do the same test on Ant/Post coordinates
  if (yLimit) {
    if (ndy > tau) {  // dissallow Posterior dir
      possibilities[LeftPostInf] = false;
      possibilities[LeftPostSup] = false;
      possibilities[RightPostInf] = false;
      possibilities[RightPostSup] = false;
    }
    else if (ndy < -tau) { // dissalow Anterior
      possibilities[LeftAntInf] = false;
      possibilities[LeftAntSup] = false;
      possibilities[RightAntInf] = false;
      possibilities[RightAntSup] = false;
    }
  }

  // Now do the same test on Inf/Sup coordinates
  if (zLimit) {
    if (ndz > tau) {  // dissallow Superior dir
      possibilities[LeftAntSup] = false;
      possibilities[LeftAntSup] = false;
      possibilities[RightAntSup] = false;
      possibilities[RightAntSup] = false;
    }
    else if (ndz < -tau) { // dissalow Inferior
      possibilities[LeftAntInf] = false;
      possibilities[LeftAntInf] = false;
      possibilities[RightAntInf] = false;
      possibilities[RightAntInf] = false;
    }
  }

  for (int child = 0; child <8; child++) {
    if (possibilities[child] && children[child]) {
      children[child]->getDistanceSqr(point, normalizedPoint, tau, distance);
    }
  }
}
